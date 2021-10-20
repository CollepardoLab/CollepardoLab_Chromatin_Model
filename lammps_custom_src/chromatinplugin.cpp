#include "lammpsplugin.h"
#include "version.h"
#include <cstring>

#include "bond_harmonic_DNA.h"
#include "hremd_steve.h"
#include "pair_aniso.h"
#include "pair_ljlambda.h"
#include "command.h"

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
  plugin.info = "Replica-exchange molecular dynamics";
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
}
