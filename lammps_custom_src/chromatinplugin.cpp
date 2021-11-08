#include "lammpsplugin.h"
#include "version.h"
#include <cstring>

#include "bond_harmonic_DNA.h"
#include "hremd_steve.h"
#include "pair_aniso.h"
#include "pair_ljlambda.h"
#include "command.h"
#include "temper_2.h"
#include "temper_2_mod.h"
#include "RIGID/fix_rigid_nve_small.h"

using namespace LAMMPS_NS;

static Bond *bondharmoniccreator(LAMMPS *lmp)
{
  return new BondHarmonic_DNA(lmp);
}

static Command *hremdstevecreator(LAMMPS *lmp)
{
  return new HremdSteve(lmp);
}

static Pair *pairanisocreator(LAMMPS *lmp)
{
  return new PairAniso(lmp);
}

static Pair *pairljlambdacreator(LAMMPS *lmp)
{
  return new PairLJLambda(lmp);
}

static Command *temper2creator(LAMMPS *lmp)
{
	return new Temper2(lmp);
}

static Command *temper2modcreator(LAMMPS *lmp)
{
	return new Temper2Mod(lmp);
}

static Fix *rigidnvesmallcreator(LAMMPS *lmp, int argc, char **argv)
{
  return new FixRigidNVESmall(lmp,argc,argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // bond harmonic style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "bond";
  plugin.name = "harmonic/DNA";
  plugin.info = "DNA rigid basepair interaction";
  plugin.author = "Steven Farr";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &bondharmoniccreator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  // REMD
  plugin.style = "command";
  plugin.name = "hremd_steve";
  plugin.info = "Hamiltonian replica-exchange molecular dynamics";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &hremdstevecreator;
  (*register_plugin)(&plugin, lmp);

  // pair aniso
  plugin.style = "pair";
  plugin.name = "aniso";
  plugin.info = "Anisotropic interaction between ellipsoids";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pairanisocreator;
  (*register_plugin)(&plugin, lmp);

  // LJ
  plugin.style = "pair";
  plugin.name = "ljlambda";
  plugin.info = "Lennard-jones interaction with switching function";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &pairljlambdacreator;
  (*register_plugin)(&plugin, lmp);

  // Temper2
  plugin.style = "command";
  plugin.name = "temper2";
  plugin.info = "Temperature replia-exchange molecular dynamics for 2 thermostats";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &temper2creator;
  (*register_plugin)(&plugin, lmp);

  // Temper2Mod
  plugin.style = "command";
  plugin.name = "temper2_mod";
  plugin.info = "Temperature replia-exchange molecular dynamics for 2 thermostats (corrected rigid body angular velocities)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &temper2modcreator;
  (*register_plugin)(&plugin, lmp);

  // FixRigidNVESmall
  //plugin.style = "fix";
  //plugin.name = "rigid/nve/small";
  //plugin.info = "Modified rigid fix";
  //plugin.author = "Trung Nguyen";
  //plugin.creator.v2 = (lammpsplugin_factory2 *) &rigidnvesmallcreator;
  //(*register_plugin)(&plugin, lmp);

}
