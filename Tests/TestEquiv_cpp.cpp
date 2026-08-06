// Equivalence tests guarding the RunCells performance work.
//
// These tests pin the floating-point and iteration-order semantics of the
// data structures on the RunCells hot path.  Each "reference_*" function is
// a verbatim copy of the production algorithm as of the start of the perf
// work; production code may be restructured for speed, but must keep
// bit-identical results and iteration order against these references.

#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "NameDouble.h"
#include "Solution.h"

namespace
{

// ----------------------------------------------------------------------
// Deterministic PRNG (no <random> to keep results identical across
// standard-library implementations).
struct LCG
{
	unsigned long long state;
	explicit LCG(unsigned long long seed) : state(seed) {}
	unsigned long long next()
	{
		state = state * 6364136223846793005ULL + 1442695040888963407ULL;
		return state >> 17;
	}
	// uniform in [0, n)
	size_t index(size_t n) { return (size_t)(next() % (unsigned long long)n); }
	// log-uniform magnitude in [1e-30, 1e0], mixed sign
	double value()
	{
		double u = (double)(next() % 1000000ULL) / 1000000.0; // [0,1)
		double mag = std::pow(10.0, -30.0 * u);
		return (next() % 5ULL == 0ULL) ? -mag : mag;
	}
};

// ----------------------------------------------------------------------
// Reference copy of cxxNameDouble::Simplify_redox (NameDouble.cxx) at the
// start of the perf work.
cxxNameDouble
reference_simplify_redox(const cxxNameDouble &nd)
{
	std::basic_string < char >::size_type indexCh;
	cxxNameDouble new_totals;
	new_totals.type = cxxNameDouble::ND_ELT_MOLES;
	{
		std::string current_ename;
		std::string const *ename_ptr;
		cxxNameDouble::const_iterator it;

		// make list of elements in new_totals
		for (it = nd.begin(); it != nd.end(); ++it)
		{
			current_ename = it->first;
			if (it->first.size() < 4)
			{
				ename_ptr = &(it->first);
			}
			else
			{
				indexCh = it->first.find("(");
				if (indexCh != std::string::npos)
				{
					current_ename = it->first.substr(0, indexCh);
					ename_ptr = &(current_ename);
				}
				else
				{
					ename_ptr = &(it->first);
				}
			}
			if (current_ename == "H" || current_ename == "O" || current_ename == "Charge")
				continue;
			new_totals[*ename_ptr] = 0;
		}
	}

	// sum totals for elements
	{
		cxxNameDouble::const_iterator old_it = nd.begin();
		cxxNameDouble::iterator new_it = new_totals.begin();
		std::string old_ename;
		std::string const *old_ename_ptr;
		while (old_it != nd.end() && new_it != new_totals.end())
		{
			if (old_it->first.size() < 4)
			{
				old_ename_ptr = &old_it->first;
			}
			else
			{
				indexCh = old_it->first.find("(");
				if (indexCh != std::string::npos)
				{
					old_ename = old_it->first.substr(0, indexCh);
					old_ename_ptr = &old_ename;
				}
				else
				{
					old_ename_ptr = &old_it->first;
				}
			}
			int j = strcmp(new_it->first.c_str(), old_ename_ptr->c_str());
			if (j < 0)
			{
				new_it++;
			}
			else if (j == 0)
			{
				new_it->second += old_it->second;
				old_it++;
			}
			else
			{
				old_it++;
			}
		}
	}
	return new_totals;
}

// ----------------------------------------------------------------------
// Reference copy of cxxSolution::Update(const cxxNameDouble &) (Solution.cxx)
// at the start of the perf work, operating on explicit maps instead of a
// cxxSolution.  Modifies master_activity in place and returns the new totals.
cxxNameDouble
reference_update(const cxxNameDouble &totals, cxxNameDouble &master_activity,
				 const cxxNameDouble &const_nd)
{
	cxxNameDouble simple_original = reference_simplify_redox(totals);
	cxxNameDouble simple_new = reference_simplify_redox(const_nd);

	cxxNameDouble factors;
	{
		// make factors
		cxxNameDouble::iterator it = simple_new.begin();
		cxxNameDouble::iterator jit = simple_original.begin();
		while (it != simple_new.end() && jit != simple_original.end())
		{
			int j = strcmp(it->first.c_str(), jit->first.c_str());
			if (j < 0)
			{
				it++;
			}
			else if (j == 0)
			{
				if (jit->second != it->second)
				{
					if (it->second > 0 && jit->second > 0)
					{
						factors[it->first] = log10(it->second / jit->second);
					}
				}
				it++;
				jit++;
			}
			else
			{
				jit++;
			}
		}
		// Now add log factors to log activities
		{
			cxxNameDouble::iterator activity_it = master_activity.begin();
			cxxNameDouble::iterator factors_it = factors.begin();
			std::string activity_ename;
			std::basic_string < char >::size_type indexCh;
			while (activity_it != master_activity.end() && factors_it != factors.end())
			{
				activity_ename = activity_it->first;
				if (factors_it->first[0] < activity_ename[0])
				{
					factors_it++;
					continue;
				}
				else if (factors_it->first[0] > activity_ename[0])
				{
					activity_it++;
					continue;
				}
				if (activity_ename.size() > 3)
				{
					indexCh = activity_ename.find("(");
					if (indexCh != std::string::npos)
					{
						activity_ename = activity_ename.substr(0, indexCh);
					}
				}
				int j = strcmp(factors_it->first.c_str(), activity_ename.c_str());
				if (j < 0)
				{
					factors_it++;
				}
				else if (j == 0)
				{
					activity_it->second += factors_it->second;
					activity_it++;
				}
				else
				{
					activity_it++;
				}
			}
		}
	}
	return simple_new;
}

// Bit-for-bit map comparison (keys, iteration order, and double bits).
bool
maps_identical(const cxxNameDouble &a, const cxxNameDouble &b)
{
	if (a.size() != b.size())
		return false;
	cxxNameDouble::const_iterator ia = a.begin(), ib = b.begin();
	for (; ia != a.end(); ++ia, ++ib)
	{
		if (ia->first != ib->first)
			return false;
		if (std::memcmp(&ia->second, &ib->second, sizeof(double)) != 0)
			return false;
	}
	return true;
}

const char *KEY_POOL[] = {
	"C(-4)", "C(4)", "Ca", "Charge", "Fe", "Fe(2)", "Fe(3)", "H", "H(0)",
	"Na", "N(3)", "N(5)", "O", "O(0)", "S(-2)", "S(6)", "U(4)", "U(6)",
	"Al", "K", "Mg", "Si", "Cl", "X",
};
const size_t KEY_POOL_N = sizeof(KEY_POOL) / sizeof(KEY_POOL[0]);

cxxNameDouble
random_name_double(LCG &rng, size_t max_keys)
{
	cxxNameDouble nd;
	size_t n = 1 + rng.index(max_keys);
	for (size_t i = 0; i < n; i++)
	{
		nd.add(KEY_POOL[rng.index(KEY_POOL_N)], rng.value());
	}
	return nd;
}

// ----------------------------------------------------------------------
// 1. Duplicate-key folding in cxxNameDouble::add is left-associative in
//    call order.
int
test_namedouble_add_order()
{
	{
		const double a = 0.1, b = 0.2, c = 0.3;
		cxxNameDouble nd;
		nd.add("Fe", a);
		nd.add("Fe", b);
		nd.add("Fe", c);
		const double expected = (a + b) + c;
		assert(nd.size() == 1);
		assert(std::memcmp(&nd.begin()->second, &expected, sizeof(double)) == 0);
	}
	{
		const double a = 1e16, b = 1.0, c = -1e16;
		cxxNameDouble nd;
		nd.add("S(6)", a);
		nd.add("S(6)", b);
		nd.add("S(6)", c);
		const double expected = (a + b) + c; // == 0.0, NOT 1.0
		assert(std::memcmp(&nd.begin()->second, &expected, sizeof(double)) == 0);
	}
	return 1;
}

// ----------------------------------------------------------------------
// 2. Production Simplify_redox is bit-identical to the reference copy on
//    randomized redox-bearing maps.
int
test_simplify_redox_reference()
{
	LCG rng(20260806ULL);
	for (int trial = 0; trial < 1000; trial++)
	{
		cxxNameDouble nd = random_name_double(rng, 16);
		cxxNameDouble prod = nd.Simplify_redox();
		cxxNameDouble ref = reference_simplify_redox(nd);
		assert(maps_identical(prod, ref));
	}
	return 1;
}

// ----------------------------------------------------------------------
// 3. Production cxxSolution::Update matches the reference algorithm on
//    synthetic solutions (totals and master_activity, bit-for-bit).
int
test_solution_update_reference()
{
	LCG rng(424242ULL);
	for (int trial = 0; trial < 200; trial++)
	{
		cxxNameDouble totals = random_name_double(rng, 12);
		cxxNameDouble new_totals = random_name_double(rng, 12);
		// master_activity carries log-activity-like values; include element
		// keys, redox keys, and short (<4 char) keys via the shared pool.
		cxxNameDouble activity = random_name_double(rng, 12);

		cxxSolution soln;
		soln.Get_totals() = totals;
		soln.Get_master_activity() = activity;

		cxxNameDouble ref_activity = activity;
		cxxNameDouble ref_totals = reference_update(totals, ref_activity, new_totals);

		soln.Update(new_totals);

		assert(maps_identical(soln.Get_totals(), ref_totals));
		assert(maps_identical(soln.Get_master_activity(), ref_activity));
	}
	return 1;
}

// ----------------------------------------------------------------------
// 4. Copy/move of cxxSolution preserves contents (dump_raw equality) and
//    self-assignment is safe.  Before the move-semantics patch, std::move
//    falls back to a copy; the assertions hold either way.
std::string
dump_solution(const cxxSolution &s)
{
	std::ostringstream oss;
	s.dump_raw(oss, 0);
	return oss.str();
}

int
test_solution_copy_move()
{
	LCG rng(77ULL);
	cxxSolution s1;
	s1.Set_n_user_both(7);
	s1.Set_tc(12.5);
	s1.Set_ph(7.25);
	s1.Set_pe(4.0);
	s1.Set_total_h(111.0);
	s1.Set_total_o(55.5);
	s1.Set_mass_water(1.0);
	s1.Get_totals() = random_name_double(rng, 12);
	s1.Get_master_activity() = random_name_double(rng, 12);

	const std::string d1 = dump_solution(s1);

	// copy-assign
	cxxSolution s2;
	s2 = s1;
	assert(dump_solution(s2) == d1);

	// copy-construct
	cxxSolution s3(s1);
	assert(dump_solution(s3) == d1);

	// move-assign (copy fallback pre-patch)
	cxxSolution s4;
	s4 = std::move(s3);
	assert(dump_solution(s4) == d1);

	// move-construct from another copy
	cxxSolution s5(s1);
	cxxSolution s6(std::move(s5));
	assert(dump_solution(s6) == d1);

	// self-assignment must be a no-op
	cxxSolution *p2 = &s2;
	s2 = *p2;
	assert(dump_solution(s2) == d1);

	// original untouched
	assert(dump_solution(s1) == d1);
	return 1;
}

// ----------------------------------------------------------------------
// 5. cxxNameDouble iterates in byte-ascending key order regardless of
//    insertion order.
int
test_namedouble_sorted_iteration()
{
	LCG rng(999ULL);
	for (int trial = 0; trial < 100; trial++)
	{
		// random insertion order over a shuffled pool subset
		cxxNameDouble nd;
		for (size_t i = 0; i < KEY_POOL_N; i++)
		{
			nd.add(KEY_POOL[rng.index(KEY_POOL_N)], rng.value());
		}
		const char *prev = NULL;
		for (cxxNameDouble::const_iterator it = nd.begin(); it != nd.end(); ++it)
		{
			if (prev != NULL)
			{
				assert(strcmp(prev, it->first.c_str()) < 0);
			}
			prev = it->first.c_str();
		}
	}
	return 1;
}

} // namespace

int
TestEquiv_cpp()
{
	test_namedouble_add_order();
	std::cerr << "TestEquiv: add order OK" << std::endl;
	test_simplify_redox_reference();
	std::cerr << "TestEquiv: Simplify_redox reference OK" << std::endl;
	test_solution_update_reference();
	std::cerr << "TestEquiv: Solution::Update reference OK" << std::endl;
	test_solution_copy_move();
	std::cerr << "TestEquiv: Solution copy/move OK" << std::endl;
	test_namedouble_sorted_iteration();
	std::cerr << "TestEquiv: sorted iteration OK" << std::endl;
	return 1;
}
