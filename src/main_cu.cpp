#include <chrono>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include "cu_helper.h"
#include "io_service_pool.hpp"
#include "safe_class.hpp"
#include "random_forest.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "log.hpp"
#include "source.hpp"

//! Represents a data wrapper for kernel callback.
template <typename T>
class callback_data
{
public:
	callback_data(io_service_pool& io, const path& output_folder_path, const size_t max_conformations, const size_t num_tasks, const receptor& rec, const forest& f, const scoring_function& sf, const T dev, const float* const cnfh, ligand&& lig_, safe_function& safe_print, log_engine& log, safe_vector<T>& idle) : io(io), output_folder_path(output_folder_path), max_conformations(max_conformations), num_tasks(num_tasks), rec(rec), f(f), sf(sf), dev(dev), cnfh(cnfh), lig(move(lig_)), safe_print(safe_print), log(log), idle(idle) {}
	io_service_pool& io;
	const path& output_folder_path;
	const size_t max_conformations;
	const size_t num_tasks;
	const receptor& rec;
	const forest& f;
	const scoring_function& sf;
	const T dev;
	const float* const cnfh;
	ligand lig;
	safe_function& safe_print;
	log_engine& log;
	safe_vector<T>& idle;
};

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path;
	array<float, 3> center, size;
	size_t seed, num_threads, num_trees, num_tasks, num_bfgs_iterations, max_conformations;
	float granularity;

	// Parse program options in a try/catch block.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_num_trees = 128;
		const size_t default_num_tasks = 256;
		const size_t default_num_bfgs_iterations = 300;
		const size_t default_max_conformations = 9;
		const  float default_granularity = 0.15625f;

		// Set up options description.
		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("receptor", value<path>(&receptor_path)->required(), "receptor in PDBQT format")
			("input_folder", value<path>(&input_folder_path)->required(), "folder of input ligands in PDBQT format")
			("center_x", value<float>(&center[0])->required(), "x coordinate of the search space center")
			("center_y", value<float>(&center[1])->required(), "y coordinate of the search space center")
			("center_z", value<float>(&center[2])->required(), "z coordinate of the search space center")
			("size_x", value<float>(&size[0])->required(), "size in the x dimension in Angstrom")
			("size_y", value<float>(&size[1])->required(), "size in the y dimension in Angstrom")
			("size_z", value<float>(&size[2])->required(), "size in the z dimension in Angstrom")
			;
		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output ligands in PDBQT format")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file in csv format")
			;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "worker threads to use")
			("trees", value<size_t>(&num_trees)->default_value(default_num_trees), "trees in random forest")
			("tasks", value<size_t>(&num_tasks)->default_value(default_num_tasks), "Monte Carlo tasks for global search")
			("generations", value<size_t>(&num_bfgs_iterations)->default_value(default_num_bfgs_iterations), "generations in BFGS")
			("max_conformations", value<size_t>(&max_conformations)->default_value(default_max_conformations), "maximum binding conformations to write")
			("granularity", value<float>(&granularity)->default_value(default_granularity), "density of probe atoms of grid maps")
			("help", "help information")
			("version", "version information")
			("config", value<path>(), "configuration file to load options from")
			;
		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);

		// If no command line argument is supplied or help is requested, print the usage and exit.
		if (argc == 1 || vm.count("help"))
		{
			cout << all_options;
			return 0;
		}

		// If version is requested, print the version and exit.
		if (vm.count("version"))
		{
			cout << "3.0.0" << endl;
			return 0;
		}

		// If a configuration file is presented, parse it.
		if (vm.count("config"))
		{
			boost::filesystem::ifstream config_file(vm["config"].as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}

		// Notify the user of parsing errors, if any.
		vm.notify();

		// Validate receptor.
		if (!is_regular_file(receptor_path))
		{
			cerr << "Receptor " << receptor_path << " does not exist or is not a regular file" << endl;
			return 1;
		}

		// Validate input_folder.
		if (!is_directory(input_folder_path))
		{
			cerr << "Input folder " << input_folder_path << " does not exist or is not a directory" << endl;
			return 1;
		}

		// Validate output_folder.
		if (exists(output_folder_path))
		{
			if (!is_directory(output_folder_path))
			{
				cerr << "Output folder " << output_folder_path << " is not a directory" << endl;
				return 1;
			}
		}
		else
		{
			if (!create_directories(output_folder_path))
			{
				cerr << "Failed to create output folder " << output_folder_path << endl;
				return 1;
			}
		}
	}
	catch (const exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;

	cout << "Creating an io service pool of " << num_threads << " worker threads" << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;
	safe_function safe_print;

	cout << "Precalculating a scoring function of " << scoring_function::n << " atom types in parallel" << endl;
	scoring_function sf;
	cnt.init((sf.n + 1) * sf.n >> 1);
	for (size_t t1 = 0; t1 < sf.n; ++t1)
	for (size_t t0 = 0; t0 <=  t1; ++t0)
	{
		io.post([&, t0, t1]()
		{
			sf.precalculate(t0, t1);
			cnt.increment();
		});
	}
	cnt.wait();

	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, granularity);

	cout << "Detecting CUDA devices with compute capability 1.1 or greater" << endl;
	checkCudaErrors(cuInit(0));
	int num_devices;
	checkCudaErrors(cuDeviceGetCount(&num_devices));
	cout << "D               Name  CC SM GMEM(MB) SMEM(KB) CMEM(KB) MAPHOST ECC TIMEOUT MODE" << endl;
	vector<CUdevice> devices;
	devices.reserve(num_devices);
	for (int dev = 0; dev < num_devices; ++dev)
	{
		// Get a device handle from an ordinal.
		CUdevice device;
		checkCudaErrors(cuDeviceGet(&device, dev));

		// Filter devices with compute capability 1.1 or greater, which is required by cuMemHostGetDevicePointer and cuStreamAddCallback.
		int major;
		int minor;
		checkCudaErrors(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device));
		checkCudaErrors(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device));
		if (major == 1 && minor == 0) continue;

		// Save the device handle.
		devices.push_back(device);

		// Get and print device attributes.
		char name[256];
		size_t totalGlobalMem;
		int multiProcessorCount;
		int sharedMemPerBlock;
		int totalConstMem;
		int canMapHostMemory;
		int ECCEnabled;
		int kernelExecTimeoutEnabled;
		int computeMode;
		checkCudaErrors(cuDeviceGetName(name, sizeof(name), device));
		checkCudaErrors(cuDeviceTotalMem(&totalGlobalMem, device));
		checkCudaErrors(cuDeviceGetAttribute(&multiProcessorCount, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, device));
		checkCudaErrors(cuDeviceGetAttribute(&sharedMemPerBlock, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, device));
		checkCudaErrors(cuDeviceGetAttribute(&totalConstMem, CU_DEVICE_ATTRIBUTE_TOTAL_CONSTANT_MEMORY, device));
		checkCudaErrors(cuDeviceGetAttribute(&canMapHostMemory, CU_DEVICE_ATTRIBUTE_CAN_MAP_HOST_MEMORY, device));
		checkCudaErrors(cuDeviceGetAttribute(&ECCEnabled, CU_DEVICE_ATTRIBUTE_ECC_ENABLED, device));
		checkCudaErrors(cuDeviceGetAttribute(&kernelExecTimeoutEnabled, CU_DEVICE_ATTRIBUTE_KERNEL_EXEC_TIMEOUT, device));
		checkCudaErrors(cuDeviceGetAttribute(&computeMode, CU_DEVICE_ATTRIBUTE_COMPUTE_MODE, device));
		cout << dev << setw(19) << name << setw(2) << major << '.' << minor << setw(3) << multiProcessorCount << setw(9) << totalGlobalMem / 1048576 << setw(9) << sharedMemPerBlock / 1024 << setw(9) << totalConstMem / 1024 << setw(8) << canMapHostMemory << setw(4) << ECCEnabled << setw(8) << kernelExecTimeoutEnabled << setw(5) << computeMode << endl;
	}
	num_devices = devices.size();
	if (!num_devices)
	{
		cerr << "No CUDA devices with compute capability 1.1 or greater detected" << endl;
		return 2;
	}

	cout << "Creating contexts and compiling kernel source for " << num_devices << " devices" << endl;
	source src;
	vector<CUcontext> contexts(num_devices);
	vector<CUfunction> functions(num_devices);
	vector<array<CUdeviceptr, sf.n>> mpsd(num_devices);
	vector<CUdeviceptr> mpsv(num_devices);
	vector<CUdeviceptr> slnv(num_devices);
	vector<CUdeviceptr> ligv(num_devices);
	vector<int*> ligh(num_devices);
	vector<CUdeviceptr> ligd(num_devices);
	vector<CUdeviceptr> slnd(num_devices);
	vector<float*> cnfh(num_devices);
	vector<size_t> lig_elems(num_devices, 2601);
	vector<size_t> sln_elems(num_devices, 3438);
	vector<size_t> cnf_elems(num_devices,   43);
	for (int dev = 0; dev < num_devices; ++dev)
	{
		// Create a context for the current device.
		checkCudaErrors(cuCtxCreate(&contexts[dev], CU_CTX_SCHED_AUTO/*CU_CTX_SCHED_YIELD*/ | CU_CTX_MAP_HOST, devices[dev]));
//		checkCudaErrors(cuCtxSetCacheConfig(CU_FUNC_CACHE_PREFER_L1));
//		checkCudaErrors(cuCtxSetSharedMemConfig(CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE));

		// Initialize just-in-time compilation options.
		const unsigned int num_jit_options = 2;
		array<CUjit_option, num_jit_options> jit_keys =
		{
			CU_JIT_MAX_REGISTERS,
			CU_JIT_CACHE_MODE
		};
		array<void*, num_jit_options> jit_vals =
		{
			(void*)32,
			(void*)CU_JIT_CACHE_OPTION_NONE // CU_JIT_CACHE_OPTION_CG, CU_JIT_CACHE_OPTION_CA
		};

		// Load the module into the current context.
		CUmodule module;
		checkCudaErrors(cuModuleLoadDataEx(&module, src.data(), num_jit_options, jit_keys.data(), jit_vals.data()));

		// Get functions from module.
		checkCudaErrors(cuModuleGetFunction(&functions[dev], module, "monte_carlo"));

		// Get symbols from module.
		CUdeviceptr sfec;
		CUdeviceptr sfdc;
		CUdeviceptr sfsc;
		CUdeviceptr cr0c;
		CUdeviceptr cr1c;
		CUdeviceptr nprc;
		CUdeviceptr gric;
		CUdeviceptr mpsc;
		CUdeviceptr nbic;
		CUdeviceptr sedc;
		CUdeviceptr slnc;
		CUdeviceptr ligc;
		size_t sfes;
		size_t sfds;
		size_t sfss;
		size_t cr0s;
		size_t cr1s;
		size_t nprs;
		size_t gris;
		size_t mpss;
		size_t nbis;
		size_t seds;
		size_t slns;
		size_t ligs;
		checkCudaErrors(cuModuleGetGlobal(&sfec, &sfes, module, "sfe")); //   8 const float*
		checkCudaErrors(cuModuleGetGlobal(&sfdc, &sfds, module, "sfd")); //   8 const float*
		checkCudaErrors(cuModuleGetGlobal(&sfsc, &sfss, module, "sfs")); //   4 int
		checkCudaErrors(cuModuleGetGlobal(&cr0c, &cr0s, module, "cr0")); //  12 float3
		checkCudaErrors(cuModuleGetGlobal(&cr1c, &cr1s, module, "cr1")); //  12 float3
		checkCudaErrors(cuModuleGetGlobal(&nprc, &nprs, module, "npr")); //  12 int3
		checkCudaErrors(cuModuleGetGlobal(&gric, &gris, module, "gri")); //   4 float
		checkCudaErrors(cuModuleGetGlobal(&mpsc, &mpss, module, "mps")); // 120 conat float* [15]
		checkCudaErrors(cuModuleGetGlobal(&nbic, &nbis, module, "nbi")); //   4 int
		checkCudaErrors(cuModuleGetGlobal(&sedc, &seds, module, "sed")); //   8 unsigned long
		checkCudaErrors(cuModuleGetGlobal(&slnc, &slns, module, "s0e")); //   8 float*
		checkCudaErrors(cuModuleGetGlobal(&ligc, &ligs, module, "lig")); //   8 const int*

		// Initialize symbols for scoring function.
		CUdeviceptr sfed;
		CUdeviceptr sfdd;
		const int sfsh = sf.ns;
		assert(sfes == sizeof(sfed));
		assert(sfds == sizeof(sfdd));
		assert(sfss == sizeof(sfsh));
		const size_t sfe_bytes = sizeof(float) * sf.e.size();
		const size_t sfd_bytes = sizeof(float) * sf.d.size();
		checkCudaErrors(cuMemAlloc(&sfed, sfe_bytes));
		checkCudaErrors(cuMemAlloc(&sfdd, sfd_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfed, sf.e.data(), sfe_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfdd, sf.d.data(), sfd_bytes));
		checkCudaErrors(cuMemcpyHtoD(sfec, &sfed, sfes));
		checkCudaErrors(cuMemcpyHtoD(sfdc, &sfdd, sfds));
		checkCudaErrors(cuMemcpyHtoD(sfsc, &sfsh, sfss));

		// Initialize symbols for receptor.
		assert(cr0s == sizeof(rec.corner0));
		assert(cr1s == sizeof(rec.corner1));
		assert(nprs == sizeof(rec.num_probes));
		assert(gris == sizeof(rec.granularity_inverse));
		assert(mpss == sizeof(float*) * sf.n);
		checkCudaErrors(cuMemcpyHtoD(cr0c, rec.corner0.data(), cr0s));
		checkCudaErrors(cuMemcpyHtoD(cr1c, rec.corner1.data(), cr1s));
		checkCudaErrors(cuMemcpyHtoD(nprc, rec.num_probes.data(), nprs));
		checkCudaErrors(cuMemcpyHtoD(gric, &rec.granularity_inverse, gris));
		mpsv[dev] = mpsc;

		// Initialize symbols for program control.
		const int nbih = num_bfgs_iterations;
		assert(nbis == sizeof(nbih));
		assert(seds == sizeof(seed));
		checkCudaErrors(cuMemcpyHtoD(nbic, &nbih, nbis));
		checkCudaErrors(cuMemcpyHtoD(sedc, &seed, seds));

		// Allocate ligh, ligd, slnd and cnfh.
		checkCudaErrors(cuMemHostAlloc((void**)&ligh[dev], sizeof(int) * lig_elems[dev], CU_MEMHOSTALLOC_DEVICEMAP));
		checkCudaErrors(cuMemHostGetDevicePointer(&ligd[dev], ligh[dev], 0));
		checkCudaErrors(cuMemAlloc(&slnd[dev], sizeof(float) * sln_elems[dev] * num_tasks));
		checkCudaErrors(cuMemHostAlloc((void**)&cnfh[dev], sizeof(float) * cnf_elems[dev] * num_tasks, 0));

		// Initialize symbols for sln and lig.
		assert(slns == sizeof(slnd[dev]));
		assert(ligs == sizeof(ligd[dev]));
		checkCudaErrors(cuMemcpyHtoD(slnc, &slnd[dev], slns));
		checkCudaErrors(cuMemcpyHtoD(ligc, &ligd[dev], ligs));
		slnv[dev] = slnc;
		ligv[dev] = ligc;

		// Pop the current context.
		checkCudaErrors(cuCtxPopCurrent(NULL));
	}
	src.clear();
	sf.clear();

	// Initialize a vector of idle devices.
	safe_vector<int> idle(num_devices);
	iota(idle.begin(), idle.end(), 0);

	cout << "Training a random forest of " << num_trees << " trees in parallel" << endl;
	forest f(num_trees, seed);
	cnt.init(num_trees);
	for (size_t i = 0; i < num_trees; ++i)
	{
		io.post([&, i]()
		{
			f[i].train(4, f.u01_s);
			cnt.increment();
		});
	}
	cnt.wait();
	f.clear();

	// Perform docking for each ligand in the input folder.
	log_engine log;
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Executing " << num_tasks << " optimization runs of " << num_bfgs_iterations << " BFGS iterations in parallel" << endl
	     << "   Index        Ligand D  pKd 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
	for (directory_iterator dir_iter(input_folder_path), const_dir_iter; dir_iter != const_dir_iter; ++dir_iter)
	{
    // Filter files with .pdbqt extension name.
		const path& input_ligand_path = dir_iter->path();
		if (input_ligand_path.extension() != ".pdbqt") continue;

    cout << "Loading ligand: " << input_ligand_path.filename() << endl;

    ligand lig(input_ligand_path);

    try {
      // Parse the ligand. Don't declare it const as it will be moved to the callback data wrapper.
      lig.load_from_path(input_ligand_path);
    } catch (const exception& e) {
      cerr << "Error loading ligand: " << input_ligand_path.filename() << endl;
      cerr << e.what() << endl;

      // Save the result with the affinities as 0
      string stem = input_ligand_path.filename().stem().string();
      lig.affinities = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      log.push_back(new log_record(move(stem), move(lig.affinities)));
      continue;
    }

		// Find atom types that are presented in the current ligand but not presented in the grid maps.
		vector<size_t> xs;
		for (size_t t = 0; t < sf.n; ++t)
		{
			if (lig.xs[t] && rec.maps[t].empty())
			{
				rec.maps[t].resize(rec.num_probes_product);
				xs.push_back(t);
			}
		}

		// Create grid maps on the fly if necessary.
		if (xs.size())
		{
			// Precalculate p_offset.
			rec.precalculate(sf, xs);

			// Create grid maps in parallel.
			cnt.init(rec.num_probes[2]);
			for (size_t z = 0; z < rec.num_probes[2]; ++z)
			{
				io.post([&,z]()
				{
					rec.populate(xs, z, sf);
					cnt.increment();
				});
			}
			cnt.wait();
		}

		// Wait until a device is ready for execution.
		const int dev = idle.safe_pop_back();

		// Push the context of the chosen device.
		checkCudaErrors(cuCtxPushCurrent(contexts[dev]));

		// Copy grid maps from host memory to device memory if necessary.
		for (size_t t = 0; t < sf.n; ++t)
		{
			if (lig.xs[t] && !mpsd[dev][t])
			{
				checkCudaErrors(cuMemAlloc(&mpsd[dev][t], rec.map_bytes));
				checkCudaErrors(cuMemcpyHtoD(mpsd[dev][t], rec.maps[t].data(), rec.map_bytes));
				checkCudaErrors(cuMemcpyHtoD(mpsv[dev] + sizeof(CUdeviceptr) * t, &mpsd[dev][t], sizeof(CUdeviceptr)));
			}
		}

		// Reallocate ligh and ligd should the current ligand elements exceed the default size.
		const size_t this_lig_elems = lig.get_lig_elems();
		if (this_lig_elems > lig_elems[dev])
		{
			checkCudaErrors(cuMemFreeHost(ligh[dev]));
			lig_elems[dev] = this_lig_elems;
			checkCudaErrors(cuMemHostAlloc((void**)&ligh[dev], sizeof(int) * lig_elems[dev], CU_MEMHOSTALLOC_DEVICEMAP));
			checkCudaErrors(cuMemHostGetDevicePointer(&ligd[dev], ligh[dev], 0));
			checkCudaErrors(cuMemcpyHtoD(ligv[dev], &ligd[dev], sizeof(ligv[dev])));
		}

		// Compute the number of shared memory bytes.
		const size_t lig_bytes = sizeof(int) * lig_elems[dev];

		// Encode the current ligand.
    try {
      lig.encode(ligh[dev]);
    } catch (const exception& e) {
      cerr << "Error encoding ligand: " << input_ligand_path.filename() << endl;
      cerr << e.what() << endl;

      // Save the result with the affinities as 0
      string stem = input_ligand_path.filename().stem().string();
      lig.affinities = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      log.push_back(new log_record(move(stem), move(lig.affinities)));
      continue;
    }

		// Reallocate slnd should the current solution elements exceed the default size.
		const size_t this_sln_elems = lig.get_sln_elems();
		if (this_sln_elems > sln_elems[dev])
		{
			checkCudaErrors(cuMemFree(slnd[dev]));
			sln_elems[dev] = this_sln_elems;
			checkCudaErrors(cuMemAlloc(&slnd[dev], sizeof(float) * sln_elems[dev] * num_tasks));
			checkCudaErrors(cuMemcpyHtoD(slnv[dev], &slnd[dev], sizeof(slnv[dev])));
		}

		// Clear the solution buffer.
		checkCudaErrors(cuMemsetD32Async(slnd[dev], 0, sln_elems[dev] * num_tasks, NULL));

		// Launch kernel.
		void* params[] = { &lig.nv, &lig.nf, &lig.na, &lig.np };
		checkCudaErrors(cuLaunchKernel(functions[dev], (num_tasks - 1) / 32 + 1, 1, 1, 32, 1, 1, lig_bytes, NULL, params, NULL));

		// Reallocate cnfh should the current conformation elements exceed the default size.
		const size_t this_cnf_elems = lig.get_cnf_elems();
		if (this_cnf_elems > cnf_elems[dev])
		{
			checkCudaErrors(cuMemFreeHost(cnfh[dev]));
			cnf_elems[dev] = this_cnf_elems;
			checkCudaErrors(cuMemHostAlloc((void**)&cnfh[dev], sizeof(float) * cnf_elems[dev] * num_tasks, 0));
		}

		// Copy conformations from device memory to host memory.
		checkCudaErrors(cuMemcpyDtoHAsync(cnfh[dev], slnd[dev], sizeof(float) * cnf_elems[dev] * num_tasks, NULL));

		// Add a callback to the compute stream.
		checkCudaErrors(cuStreamAddCallback(NULL, [](CUstream stream, CUresult error, void* data)
		{
			checkCudaErrors(error);
			const shared_ptr<callback_data<int>> cbd(reinterpret_cast<callback_data<int>*>(data));
			cbd->io.post([=]()
			{
				const auto& output_folder_path = cbd->output_folder_path;
				const auto  max_conformations = cbd->max_conformations;
				const auto  num_tasks = cbd->num_tasks;
				const auto& rec = cbd->rec;
				const auto& f = cbd->f;
				const auto& sf = cbd->sf;
				const auto  dev = cbd->dev;
				const auto cnfh = cbd->cnfh;
				auto& lig = cbd->lig;
				auto& safe_print = cbd->safe_print;
				auto& log = cbd->log;
				auto& idle = cbd->idle;

				// Write conformations.
				lig.write(cnfh, output_folder_path, max_conformations, num_tasks, rec, f, sf);

				// Output and save ligand stem and predicted affinities.
				safe_print([&]()
				{
					string stem = lig.filename.stem().string();
					cout << setw(8) << log.size() + 1 << setw(14) << stem << setw(2) << dev << ' ';
					for_each(lig.affinities.cbegin(), lig.affinities.cbegin() + min<size_t>(lig.affinities.size(), 9), [](const float a)
					{
						cout << setw(6) << a;
					});
					cout << endl;
					log.push_back(new log_record(move(stem), move(lig.affinities)));
				});

				// Signal the main thread to post another task.
				idle.safe_push_back(dev);
			});
		}, new callback_data<int>(io, output_folder_path, max_conformations, num_tasks, rec, f, sf, dev, cnfh[dev], move(lig), safe_print, log, idle), 0));

		// Pop the context after use.
		checkCudaErrors(cuCtxPopCurrent(NULL));
	}

	// Synchronize contexts.
	for (auto& context : contexts)
	{
		checkCudaErrors(cuCtxPushCurrent(context));
		checkCudaErrors(cuCtxSynchronize());
		checkCudaErrors(cuCtxPopCurrent(NULL));
	}

	// Wait until the io service pool has finished all its tasks.
	io.wait();
	assert(idle.size() == num_devices);

	// Destroy contexts.
	for (auto& context : contexts)
	{
		checkCudaErrors(cuCtxDestroy(context));
	}

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
