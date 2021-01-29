/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   DNA rigid basepair interaction, modified from bond_harmonic.h
   S.FARR 9/2018
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   DNA rigid basepair interaction, modified from bond_harmonic.cpp
   S.FARR 2018
------------------------------------------------------------------------- */


#ifdef BOND_CLASS

BondStyle(harmonic/DNA,BondHarmonic_DNA)

#else

#ifndef LMP_BOND_HARMONIC_DNA_H
#define LMP_BOND_HARMONIC_DNA_H

#include <cstdio>
#include "bond.h"
#include "atom_vec_ellipsoid.h"
#include <unordered_map>
#include <vector>
#include <string>

// Conversion factor for helical paramters
// They are in Kcal/mol/Angstrom which agree with lammps "real" units
#define E_CONV_FACTOR_STEVE 1.0

namespace LAMMPS_NS {

class BondHarmonic_DNA : public Bond {
 public:
  BondHarmonic_DNA(class LAMMPS *);
  virtual ~BondHarmonic_DNA();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  virtual void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);
  virtual void *extract(char *, int &);
  void init_style();

protected:
  class AtomVecEllipsoid *avec;

  // The equilibrium values and stiffness matrices for each base-pair step
  // combination.
  struct helical_paramsstruct {
      double means[6];
      double K[6][6];

      //constructor and initialisation list
      helical_paramsstruct() : means{}, K{} {
        // required (empty) body of constructor
      }
  };

  // map linking bp step type to helical params
  std::unordered_map<std::string, helical_paramsstruct> helical_params_map;

  // map linking global atom tag of DNA to its base_pairs
  std::unordered_map<tagint,std::string> base_pairs;

  void get_helical_param_map(std::string fname);
  void get_basepairs(std::string fname);

  double compute_bond_energy(const double * x1,  double * q1, const double * x2, double * q2,struct  helical_paramsstruct * helical_params);
  double *k,*r0;

  virtual void allocate();

  // math functions

  double arcos(double x);
  double vec_dot_mat_dot_vec(const double * vec, double mat[6][6]);
  void rotation(const double *x, double *y, double *axis, double angle);
  void compute_helical_parameters(const double * x1, const double * ex1_in, const double * ey1_in, const double * ez1_in, const double * x2, const double * ex2_in, const double * ey2_in,const  double * ez2_in, double * out,double * mstx_out, double * msty_out, double * mstz_out);
  double mag_vec(const double * v);

  };

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
