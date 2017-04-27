#include <chrono>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include "io_service_pool.hpp"
#include "safe_class.hpp"
#include "random_forest.hpp"
#include "receptor.hpp"
#include "ligand.hpp"
#include "log.hpp"
#include "kernel.hpp"
using namespace boost::filesystem;

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

  boost::filesystem::ofstream err_ofs(err_path);

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;
	mt19937_64 rng(seed);

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
	sf.clear();

	cout << "Parsing receptor " << receptor_path << endl;
	receptor rec(receptor_path, center, size, granularity);

	vector<int>   ligh(2601);
	vector<float> slnd(3438 * num_tasks);

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
	     << "   Index        Ligand    pKd 1     2     3     4     5     6     7     8     9" << endl << setprecision(2);
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

		// Reallocate ligh and ligd should the current ligand elements exceed the default size.
		const size_t this_lig_elems = lig.get_lig_elems();
		if (this_lig_elems > ligh.size())
		{
			ligh.resize(this_lig_elems);
		}

		// Encode the current ligand.
		lig.encode(ligh.data());

		// Reallocate slnd should the current solution elements exceed the default size.
		const size_t this_sln_elems = lig.get_sln_elems() * num_tasks;
		if (this_sln_elems > slnd.size())
		{
			slnd.resize(this_sln_elems);
		}

		// Clear the solution buffer.
		slnd.assign(slnd.size(), 0);

		// Launch kernel.
		cnt.init(num_tasks);
		for (int gid = 0; gid < num_tasks; ++gid)
		{
			const size_t s = rng();
			io.post([&, s, gid]()
			{
				monte_carlo(slnd.data(), ligh.data(), lig.nv, lig.nf, lig.na, lig.np, s, num_bfgs_iterations, sf.e.data(), sf.d.data(), sf.ns, rec.corner0, rec.corner1, rec.num_probes, rec.granularity_inverse, rec.maps, gid, num_tasks);
				cnt.increment();
			});
		}
		cnt.wait();

		// Reallocate cnfh should the current conformation elements exceed the default size.
		const size_t this_cnf_elems = lig.get_cnf_elems() * num_tasks;

		io.post(bind([&](ligand lig, vector<float> cnfh)
		{
			// Write conformations.
			lig.write(cnfh.data(), output_folder_path, max_conformations, num_tasks, rec, f, sf);

			// Output and save ligand stem and predicted affinities.
			safe_print([&]()
			{
				string stem = lig.filename.stem().string();
				cout << setw(8) << log.size() + 1 << setw(14) << stem << setw(2) << "   ";
				for_each(lig.affinities.cbegin(), lig.affinities.cbegin() + min<size_t>(lig.affinities.size(), 9), [](const float a)
				{
					cout << setw(6) << a;
				});
				cout << endl;
				log.push_back(new log_record(move(stem), move(lig.affinities)));
			});
		}, move(lig), vector<float>(slnd.cbegin(), slnd.cbegin() + this_cnf_elems)));
	}

	// Wait until the io service pool has finished all its tasks.
	io.wait();

	// Sort and write ligand log records to the log file.
	if (log.empty()) return 0;
	cout << "Writing log records of " << log.size() << " ligands to " << log_path << endl;
	log.sort();
	log.write(log_path);
}
