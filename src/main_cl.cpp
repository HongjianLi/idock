#include <chrono>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <CL/cl.h>
#include "cl_helper.h"
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
	callback_data(io_service_pool& io, cl_event cbex, const path& output_folder_path, const size_t max_conformations, const size_t num_tasks, const receptor& rec, const forest& f, const scoring_function& sf, const T dev, float* const cnfh, ligand&& lig_, cl_mem slnd, safe_function& safe_print, log_engine& log, safe_vector<T>& idle) : io(io), cbex(cbex), output_folder_path(output_folder_path), max_conformations(max_conformations), num_tasks(num_tasks), rec(rec), f(f), sf(sf), dev(dev), cnfh(cnfh), lig(move(lig_)), slnd(slnd), safe_print(safe_print), log(log), idle(idle) {}
	io_service_pool& io;
	cl_event cbex;
	const path& output_folder_path;
	const size_t max_conformations;
	const size_t num_tasks;
	const receptor& rec;
	const forest& f;
	const scoring_function& sf;
	const T dev;
	float* const cnfh;
	ligand lig;
	cl_mem slnd;
	safe_function& safe_print;
	log_engine& log;
	safe_vector<T>& idle;
};

int main(int argc, char* argv[])
{
	path receptor_path, input_folder_path, output_folder_path, log_path, err_path;
	array<float, 3> center, size;
	size_t seed, num_threads, num_trees, num_tasks, num_bfgs_iterations, max_conformations;
	float granularity;

	// Parse program options in a try/catch block.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.csv";
    const path default_err_path = "err.csv";
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
      ("err", value<path>(&err_path)->default_value(default_err_path), "err file in csv format")
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
	const int sfs = sf.ns;

	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, granularity);

	cout << "Detecting OpenCL platforms" << endl;
	char name[256];
	char version[256];
	cl_uint num_platforms;
	checkOclErrors(clGetPlatformIDs(0, NULL, &num_platforms));
	vector<cl_platform_id> platforms(num_platforms);
	checkOclErrors(clGetPlatformIDs(num_platforms, platforms.data(), NULL));
	vector<cl_uint> num_platform_devices(num_platforms);
	cl_uint num_devices = 0;
	for (cl_uint i = 0; i < num_platforms; ++i)
	{
		const auto platform = platforms[i];
		checkOclErrors(clGetPlatformInfo(platform, CL_PLATFORM_NAME, sizeof(name), name, NULL));
		checkOclErrors(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, sizeof(version), version, NULL));
		cout << i << ' ' << name << ", " << version << endl;
		checkOclErrors(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &num_platform_devices[i]));
		num_devices += num_platform_devices[i];
	}

	cout << "Detecting OpenCL devices" << endl;
	if (!num_devices)
	{
		cerr << "No OpenCL devices detected" << endl;
		return 2;
	}
	vector<cl_device_id> devices(num_devices);
	vector<bool> cl12(num_devices);
	vector<cl_bool> host_unified_memory(num_devices);
	cout << "D                           Name  CL CU GMEM(MB) LMEM(KB) CMEM(KB) LMEMTYPE ECC" << endl;
	const char* local_mem_types[] = { "NONE", "LOCAL", "GLOBAL" };
	for (cl_uint i = 0, dev = 0; i < num_platforms; ++i)
	{
		const auto platform = platforms[i];
		const auto npd = num_platform_devices[i];
		checkOclErrors(clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, npd, &devices[dev], NULL));
		for (int d = 0; d < npd; ++d, ++dev)
		{
			const auto device = devices[dev];
			cl_uint max_compute_units;
			cl_ulong global_mem_size;
			cl_ulong local_mem_size;
			cl_ulong max_constant_buffer_size;
			cl_bool error_correction_support;
			cl_device_local_mem_type local_mem_type;
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(name), name, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION, sizeof(version), version, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(max_compute_units), &max_compute_units, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(global_mem_size), &global_mem_size, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_mem_size), &local_mem_size, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(max_constant_buffer_size), &max_constant_buffer_size, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_HOST_UNIFIED_MEMORY, sizeof(host_unified_memory[dev]), &host_unified_memory[dev], NULL));
			checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL));
			cl12[dev] = version[9] > '1' || version[11] >= '2';
			const string name_str(name);
			const size_t name_len = name_str.size();
			cout << dev << setw(31) << (name_len <= 30 ? name_str : name_str.substr(name_str.find(' ', name_len - 31) + 1)) << ' ' << version[9] << '.' << version[11] << setw(3) << max_compute_units << setw(9) << global_mem_size / 1048576 << setw(9) << local_mem_size / 1024 << setw(9) << max_constant_buffer_size / 1024 << setw(9) << local_mem_types[local_mem_type] << setw(4) << error_correction_support << endl;
		}
	}

	cout << "Creating contexts and compiling kernel source for " << num_devices << " devices" << endl;
	source src;
	const char* sources[] = { src.data() };
	const size_t source_length = src.size();
	vector<cl_context> contexts(num_devices);
	vector<cl_command_queue> queues(num_devices);
	vector<cl_program> programs(num_devices);
	vector<cl_kernel> kernels(num_devices);
	vector<cl_mem> sfed(num_devices);
	vector<cl_mem> sfdd(num_devices);
	vector<cl_mem> ligd(num_devices);
	vector<cl_mem> slnd(num_devices);
	vector<size_t> lig_elems(num_devices, 2601);
	vector<size_t> sln_elems(num_devices, 3438);
	vector<size_t> cnf_elems(num_devices,   43);
	vector<array<cl_mem, sf.n>> mpsd(num_devices);
	cl_int error;
	for (int dev = 0; dev < num_devices; ++dev)
	{
		// Get device.
		cl_device_id device = devices[dev];

		// Create context.
		cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, &error);
		checkOclErrors(error);
		contexts[dev] = context;

		// Create command queue.
		cl_command_queue_properties queue_properties;
		checkOclErrors(clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(queue_properties), &queue_properties, NULL));
		cl_command_queue queue = clCreateCommandQueue(context, device, queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE ? CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE : 0/* | CL_QUEUE_PROFILING_ENABLE*/, &error);
		checkOclErrors(error);
		queues[dev] = queue;

		// Create program.
		cl_program program = clCreateProgramWithSource(context, 1, sources, &source_length, &error);
		checkOclErrors(error);
		programs[dev] = program;

		// Build program.
		checkOclErrors(clBuildProgram(program, 0, NULL, "-cl-fast-relaxed-math"/*-cl-std=CL1.2 -cl-nv-maxrregcount 32*/, NULL, NULL));

		// Create kernel from program.
		cl_kernel kernel = clCreateKernel(program, "monte_carlo", &error);
		checkOclErrors(error);
		kernels[dev] = kernel;

		// Create buffers for sfe and sfd.
		const size_t sfe_bytes = sizeof(float) * sf.e.size();
		const size_t sfd_bytes = sizeof(float) * sf.d.size();
		sfed[dev] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sfe_bytes, sf.e.data(), &error);
		checkOclErrors(error);
		sfdd[dev] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sfd_bytes, sf.d.data(), &error);
		checkOclErrors(error);

		// Create buffers for ligh, ligd, slnd and cnfh.
		ligd[dev] = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int) * lig_elems[dev], NULL, &error);
		checkOclErrors(error);
		slnd[dev] = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * sln_elems[dev] * num_tasks, NULL, &error);
		checkOclErrors(error);

		// Set kernel arguments.
		checkOclErrors(clSetKernelArg(kernel,  0, sizeof(cl_mem), &slnd[dev]));
		checkOclErrors(clSetKernelArg(kernel,  1, sizeof(cl_mem), &ligd[dev]));
		checkOclErrors(clSetKernelArg(kernel,  7, sizeof(int), &num_bfgs_iterations));
		checkOclErrors(clSetKernelArg(kernel,  8, sizeof(cl_mem), &sfed[dev]));
		checkOclErrors(clSetKernelArg(kernel,  9, sizeof(cl_mem), &sfdd[dev]));
		checkOclErrors(clSetKernelArg(kernel, 10, sizeof(int), &sfs));
		checkOclErrors(clSetKernelArg(kernel, 11, sizeof(cl_float3), rec.corner0.data()));
		checkOclErrors(clSetKernelArg(kernel, 12, sizeof(cl_float3), rec.corner1.data()));
		checkOclErrors(clSetKernelArg(kernel, 13, sizeof(cl_float3), rec.num_probes.data()));
		checkOclErrors(clSetKernelArg(kernel, 14, sizeof(rec.granularity_inverse), &rec.granularity_inverse));
		for (size_t t = 0; t < sf.n; ++t)
		{
			checkOclErrors(clSetKernelArg(kernel, 15 + t, sizeof(cl_mem), &mpsd[dev][t]));
		}
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
	vector<cl_event> cbex(num_devices);
	cout.setf(ios::fixed, ios::floatfield);
	cout << "Executing " << num_tasks << " optimization runs of " << num_bfgs_iterations << " BFGS iterations in parallel" << endl
	     << "   Index        Ligand D  pKd 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
	for (directory_iterator dir_iter(input_folder_path), const_dir_iter; dir_iter != const_dir_iter; ++dir_iter)
	{
		// Filter files with .pdbqt extension name.
		const path& input_ligand_path = dir_iter->path();
		if (input_ligand_path.extension() != ".pdbqt") continue;

		// Parse the ligand. Don't declare it const as it will be moved to the callback data wrapper.
		ligand lig(input_ligand_path);

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

		// Copy grid maps from host memory to device memory if necessary.
		for (size_t t = 0; t < sf.n; ++t)
		{
			if (lig.xs[t] && !mpsd[dev][t])
			{
				mpsd[dev][t] = clCreateBuffer(contexts[dev], CL_MEM_READ_ONLY, rec.map_bytes, NULL, &error);
				checkOclErrors(error);
				checkOclErrors(clEnqueueWriteBuffer(queues[dev], mpsd[dev][t], CL_TRUE, 0, rec.map_bytes, rec.maps[t].data(), 0, NULL, NULL));
				checkOclErrors(clSetKernelArg(kernels[dev], 15 + t, sizeof(cl_mem), &mpsd[dev][t]));
			}
		}

		// Reallocate ligd should the current ligand elements exceed the default size.
		const size_t this_lig_elems = lig.get_lig_elems();
		if (this_lig_elems > lig_elems[dev])
		{
			checkOclErrors(clReleaseMemObject(ligd[dev]));
			lig_elems[dev] = this_lig_elems;
			ligd[dev] = clCreateBuffer(contexts[dev], CL_MEM_READ_ONLY, sizeof(int) * lig_elems[dev], NULL, &error);
			checkOclErrors(error);
			checkOclErrors(clSetKernelArg(kernels[dev], 1, sizeof(cl_mem), &ligd[dev]));
		}

		// Compute the number of local memory bytes.
		const size_t lig_bytes = sizeof(int) * lig_elems[dev];

		// Encode the current ligand.
		cl_event input_events[2];
		int* ligh = (int*)clEnqueueMapBuffer(queues[dev], ligd[dev], CL_TRUE, cl12[dev] ? CL_MAP_WRITE_INVALIDATE_REGION : CL_MAP_WRITE, 0, lig_bytes, 0, NULL, NULL, &error);
		checkOclErrors(error);
		lig.encode(ligh);
		checkOclErrors(clEnqueueUnmapMemObject(queues[dev], ligd[dev], ligh, 0, NULL, &input_events[0]));

		// Reallocate slnd should the current solution elements exceed the default size.
		const size_t this_sln_elems = lig.get_sln_elems();
		if (this_sln_elems > sln_elems[dev])
		{
			checkOclErrors(clReleaseMemObject(slnd[dev]));
			sln_elems[dev] = this_sln_elems;
			slnd[dev] = clCreateBuffer(contexts[dev], CL_MEM_READ_WRITE, sizeof(float) * sln_elems[dev] * num_tasks, NULL, &error);
			checkOclErrors(error);
			checkOclErrors(clSetKernelArg(kernels[dev], 0, sizeof(cl_mem), &slnd[dev]));
		}

		// Clear the solution buffer.
		if (cl12[dev])
		{
			const float pattern = 0.0f;
			checkOclErrors(clEnqueueFillBuffer(queues[dev], slnd[dev], &pattern, sizeof(pattern), 0, sizeof(float) * sln_elems[dev] * num_tasks, 0, NULL, &input_events[1]));
		}
		else
		{
			float* slnh = (float*)clEnqueueMapBuffer(queues[dev], slnd[dev], CL_TRUE, CL_MAP_WRITE, 0, sizeof(float) * sln_elems[dev] * num_tasks, 0, NULL, NULL, &error);
			checkOclErrors(error);
			memset(slnh, 0, sizeof(float) * sln_elems[dev] * num_tasks);
			checkOclErrors(clEnqueueUnmapMemObject(queues[dev], slnd[dev], slnh, 0, NULL, &input_events[1]));
		}

		// Launch kernel.
		checkOclErrors(clSetKernelArg(kernels[dev],  2, sizeof(int), &lig.nv));
		checkOclErrors(clSetKernelArg(kernels[dev],  3, sizeof(int), &lig.nf));
		checkOclErrors(clSetKernelArg(kernels[dev],  4, sizeof(int), &lig.na));
		checkOclErrors(clSetKernelArg(kernels[dev],  5, sizeof(int), &lig.np));
		checkOclErrors(clSetKernelArg(kernels[dev],  6, lig_bytes, NULL));
		const size_t gws = num_tasks;
		const size_t lws = 32;
		cl_event kernel_event;
		checkOclErrors(clEnqueueNDRangeKernel(queues[dev], kernels[dev], 1, NULL, &gws, &lws, 2, input_events, &kernel_event));

		// Reallocate cnfh should the current conformation elements exceed the default size.
		const size_t this_cnf_elems = lig.get_cnf_elems();
		if (this_cnf_elems > cnf_elems[dev])
		{
//			checkOclErrors(clReleaseMemObject(cnfh[dev]));
			cnf_elems[dev] = this_cnf_elems;
//			cnfh[dev] = clCreateBuffer(contexts[dev], CL_MEM_ALLOC_HOST_PTR, sizeof(float) * cnf_elems[dev] * num_tasks, NULL, &error);
//			checkOclErrors(error);
		}

		// Copy conformations from device memory to host memory.
		cl_event output_event;
		float* cnfh = (float*)clEnqueueMapBuffer(queues[dev], slnd[dev], CL_FALSE, CL_MAP_READ, 0, sizeof(float) * cnf_elems[dev] * num_tasks, 1, &kernel_event, &output_event, &error);
		checkOclErrors(error);

		// Create callback events.
		if (cbex[dev]) checkOclErrors(clReleaseEvent(cbex[dev]));
		cbex[dev] = clCreateUserEvent(contexts[dev], &error);
		checkOclErrors(error);

		// Add a callback to the output event.
		checkOclErrors(clSetEventCallback(output_event, CL_COMPLETE, [](cl_event event, cl_int command_exec_status, void* data)
		{
			assert(command_exec_status == CL_COMPLETE);
			cl_command_queue queue;
			checkOclErrors(clGetEventInfo(event, CL_EVENT_COMMAND_QUEUE, sizeof(queue), &queue, 0));
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
				auto slnd = cbd->slnd;
				auto& safe_print = cbd->safe_print;
				auto& log = cbd->log;
				auto& idle = cbd->idle;

				// Write conformations.
				lig.write(cnfh, output_folder_path, max_conformations, num_tasks, rec, f, sf);

				// Unmap cnfh.
				checkOclErrors(clEnqueueUnmapMemObject(queue, slnd, cnfh, 0, NULL, NULL));

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
			checkOclErrors(clSetUserEventStatus(cbd->cbex, CL_COMPLETE));
		}, new callback_data<int>(io, cbex[dev], output_folder_path, max_conformations, num_tasks, rec, f, sf, dev, cnfh, move(lig), slnd[dev], safe_print, log, idle)));
	}

	// Synchronize queues and callback events.
	for (int dev = 0; dev < num_devices; ++dev)
	{
		checkOclErrors(clFinish(queues[dev]));
		if (cbex[dev]) checkOclErrors(clWaitForEvents(1, &cbex[dev]));
	}

	// Wait until the io service pool has finished all its tasks.
	io.wait();
	assert(idle.size() == num_devices);

	// Release resources.
	for (int dev = 0; dev < num_devices; ++dev)
	{
		for (auto mapd : mpsd[dev])
		{
			if (mapd) checkOclErrors(clReleaseMemObject(mapd));
		}
		if (cbex[dev]) checkOclErrors(clReleaseEvent(cbex[dev]));
		checkOclErrors(clReleaseMemObject(sfdd[dev]));
		checkOclErrors(clReleaseMemObject(sfed[dev]));
		checkOclErrors(clReleaseMemObject(slnd[dev]));
		checkOclErrors(clReleaseMemObject(ligd[dev]));
		checkOclErrors(clReleaseKernel(kernels[dev]));
		checkOclErrors(clReleaseProgram(programs[dev]));
		checkOclErrors(clReleaseCommandQueue(queues[dev]));
		checkOclErrors(clReleaseContext(contexts[dev]));
	}

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
