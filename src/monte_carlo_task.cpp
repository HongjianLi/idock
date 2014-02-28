#include <random>
#include "monte_carlo_task.hpp"
using namespace std;
//using std::random::mt19937_64;

void monte_carlo_task(ptr_vector<result>& results, const ligand& lig, const size_t seed, const scoring_function& sf, const box& b, const vector<vector<double>>& grid_maps)
{
	// Define constants.
	const size_t num_mc_iterations = 100 * lig.num_heavy_atoms; ///< The number of iterations correlates to the complexity of ligand.
	const size_t num_entities  = 2 + lig.num_active_torsions; // Number of entities to mutate.
	const size_t num_variables = 6 + lig.num_active_torsions; // Number of variables to optimize.
	const double e_upper_bound = static_cast<double>(4 * lig.num_heavy_atoms); // A conformation will be droped if its free energy is not better than e_upper_bound.
	const double required_square_error = static_cast<double>(1 * lig.num_heavy_atoms); // Ligands with RMSD < 1.0 will be clustered into the same cluster.
	const double pi = static_cast<double>(3.1415926535897932); ///< Pi.

	mt19937_64 eng(seed);
	uniform_real_distribution<double> uniform_01_gen(  0,  1);
	uniform_real_distribution<double> uniform_11_gen( -1,  1);
	uniform_real_distribution<double> uniform_pi_gen(-pi, pi);
	uniform_real_distribution<double> uniform_box0_gen(b.corner1[0], b.corner2[0]);
	uniform_real_distribution<double> uniform_box1_gen(b.corner1[1], b.corner2[1]);
	uniform_real_distribution<double> uniform_box2_gen(b.corner1[2], b.corner2[2]);
	uniform_int_distribution<size_t> uniform_entity_gen(0, num_entities - 1);
	normal_distribution<double> normal_01_gen(0, 1);

	// Generate an initial random conformation c0, and evaluate it.
	conformation c0(lig.num_active_torsions);
	double e0, f0;
	change g0(lig.num_active_torsions);
	bool valid_conformation = false;
	for (size_t i = 0; (i < 1000) && (!valid_conformation); ++i)
	{
		// Randomize conformation c0.
		c0.position = array<double, 3>{uniform_box0_gen(eng), uniform_box1_gen(eng), uniform_box2_gen(eng)};
		c0.orientation = normalize(array<double, 4>{normal_01_gen(eng), normal_01_gen(eng), normal_01_gen(eng), normal_01_gen(eng)});
		for (size_t i = 0; i < lig.num_active_torsions; ++i)
		{
			c0.torsions[i] = uniform_pi_gen(eng);
		}
		valid_conformation = lig.evaluate(c0, sf, b, grid_maps, e_upper_bound, e0, f0, g0);
	}
	if (!valid_conformation) return;
	double best_e = e0; // The best free energy so far.

	// Initialize necessary variables for BFGS.
	conformation c1(lig.num_active_torsions), c2(lig.num_active_torsions); // c2 = c1 + ap.
	double e1, f1, e2, f2;
	change g1(lig.num_active_torsions), g2(lig.num_active_torsions);
	change p(lig.num_active_torsions); // Descent direction.
	double alpha, pg1, pg2; // pg1 = p * g1. pg2 = p * g2.
	size_t num_alpha_trials;

	// Initialize the inverse Hessian matrix to identity matrix.
	// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
	// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
	// See N&R for a recipe to find this initializer.
	triangular_matrix<double> identity_hessian(num_variables, 0); // Symmetric triangular matrix.
	for (size_t i = 0; i < num_variables; ++i)
		identity_hessian[triangular_matrix_restrictive_index(i, i)] = 1;

	// Initialize necessary variables for updating the Hessian matrix h.
	triangular_matrix<double> h(identity_hessian);
	change y(lig.num_active_torsions); // y = g2 - g1.
	change mhy(lig.num_active_torsions); // mhy = -h * y.
	double yhy, yp, ryp, pco;

	for (size_t mc_i = 0; mc_i < num_mc_iterations; ++mc_i)
	{
		size_t num_mutations = 0;
		size_t mutation_entity;

		// Mutate c0 into c1, and evaluate c1.
		do
		{
			// Make a copy, so the previous conformation is retained.
			c1 = c0;

			// Determine an entity to mutate.
			mutation_entity = uniform_entity_gen(eng);
			assert(mutation_entity < num_entities);
			if (mutation_entity < lig.num_active_torsions) // Mutate an active torsion.
			{
				c1.torsions[mutation_entity] = uniform_pi_gen(eng);
			}
			else if (mutation_entity == lig.num_active_torsions) // Mutate position.
			{
				c1.position += array<double, 3>{uniform_11_gen(eng), uniform_11_gen(eng), uniform_11_gen(eng)};
			}
			else // Mutate orientation.
			{
				c1.orientation = vec3_to_qtn4(static_cast<double>(0.01) * array<double, 3>{uniform_11_gen(eng), uniform_11_gen(eng), uniform_11_gen(eng)}) * c1.orientation;
				assert(normalized(c1.orientation));
			}
			++num_mutations;
		} while (!lig.evaluate(c1, sf, b, grid_maps, e_upper_bound, e1, f1, g1));

		// Initialize the Hessian matrix to identity.
		h = identity_hessian;

		// Given the mutated conformation c1, use BFGS to find a local minimum.
		// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when an appropriate alpha cannot be found.
		while (true)
		{
			// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
			for (size_t i = 0; i < num_variables; ++i)
			{
				double sum = 0;
				for (size_t j = 0; j < num_variables; ++j)
					sum += h[triangular_matrix_permissive_index(i, j)] * g1[j];
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (size_t i = 0; i < num_variables; ++i)
				pg1 += p[i] * g1[i];

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for num_alphas times.
			// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
			alpha = 1.0;
			for (num_alpha_trials = 0; num_alpha_trials < num_alphas; ++num_alpha_trials)
			{
				// Obtain alpha from the precalculated alpha values.
				alpha *= 0.1;

				// Calculate c2 = c1 + ap.
				c2.position = c1.position + alpha * array<double, 3>{p[0], p[1], p[2]};
				assert(normalized(c1.orientation));
				c2.orientation = vec3_to_qtn4(alpha * array<double, 3>{p[3], p[4], p[5]}) * c1.orientation;
				assert(normalized(c2.orientation));
				for (size_t i = 0; i < lig.num_active_torsions; ++i)
				{
					c2.torsions[i] = c1.torsions[i] + alpha * p[6 + i];
				}

				// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (lig.evaluate(c2, sf, b, grid_maps, e1 + 0.0001 * alpha * pg1, e2, f2, g2))
				{
					pg2 = 0;
					for (size_t i = 0; i < num_variables; ++i)
						pg2 += p[i] * g2[i];
					if (pg2 >= 0.9 * pg1)
						break; // An appropriate alpha is found.
				}
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (num_alpha_trials == num_alphas) break;

			// Update Hessian matrix h.
			for (size_t i = 0; i < num_variables; ++i) // Calculate y = g2 - g1.
				y[i] = g2[i] - g1[i];
			for (size_t i = 0; i < num_variables; ++i) // Calculate mhy = -h * y.
			{
				double sum = 0;
				for (size_t j = 0; j < num_variables; ++j)
					sum += h[triangular_matrix_permissive_index(i, j)] * y[j];
				mhy[i] = -sum;
			}
			yhy = 0;
			for (size_t i = 0; i < num_variables; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
				yhy -= y[i] * mhy[i];
			yp = 0;
			for (size_t i = 0; i < num_variables; ++i) // Calculate yp = y * p.
				yp += y[i] * p[i];
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (size_t i = 0; i < num_variables; ++i)
			for (size_t j = i; j < num_variables; ++j) // includes i
			{
				h[triangular_matrix_restrictive_index(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
			}

			// Move to the next iteration.
			c1 = c2;
			e1 = e2;
			f1 = f2;
			g1 = g2;
		}

		// Accept c1 according to Metropolis criteria.
		const double delta = e0 - e1;
		if ((delta > 0) || (uniform_01_gen(eng) < exp(delta)))
		{
			// best_e is the best energy of all the conformations in the container.
			// e1 will be saved if and only if it is even better than the best one.
			if (e1 < best_e || results.size() < results.capacity())
			{
				add_to_result_container(results, lig.compose_result(e1, f1, c1), required_square_error);
				if (e1 < best_e) best_e = e0;
			}

			// Save c1 into c0.
			c0 = c1;
			e0 = e1;
		}
	}
}
