/*****************/
// Debug dump utilities for force parity checks
/*****************/
#include "force_debug.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace force_debug {
namespace {

bool enabled() {
  const char *env = std::getenv("CELLDYN_FORCE_DUMP");
  if (!env) return false;
  return std::string(env) == "1";
}

std::string out_dir() {
  const char *env = std::getenv("CELLDYN_FORCE_DUMP_DIR");
  if (!env || std::string(env).empty()) return "force_dump";
  return std::string(env);
}

bool should_dump_step(unsigned int step) {
  const char *env = std::getenv("CELLDYN_FORCE_DUMP_EVERY");
  if (!env) return true;
  int every = std::atoi(env);
  if (every <= 0) return true;
  return step % static_cast<unsigned int>(every) == 0;
}

std::string zero_pad(unsigned int value, int width = 10) {
  std::ostringstream oss;
  oss << std::setw(width) << std::setfill('0') << value;
  return oss.str();
}

} // namespace

void dump_force_snapshot(const Global *p_g, int deg, const char *phase) {
  if (!enabled()) return;
  if (!should_dump_step(p_g->step)) return;

  std::string dir = out_dir();
  std::string mk = "mkdir -p \"" + dir + "\"";
  (void)std::system(mk.c_str());

  std::string tag = phase ? std::string(phase) : std::string("phase");
  std::string f_force = dir + "/force_" + tag + "_" + zero_pad(p_g->step) + ".csv";
  std::string f_lt = dir + "/lt_" + tag + "_" + zero_pad(p_g->step) + ".csv";

  {
    std::ofstream ofs(f_force);
    ofs << "step,deg,vertex_id,fx,fy,fz\n";
    ofs << std::setprecision(17);
    for (int i = 0; i < (int)p_g->p_v.size(); i++) {
      const Vertex *vp = p_g->p_v[i];
      ofs << p_g->step << "," << deg << "," << i << ","
          << vp->frc[deg].x << "," << vp->frc[deg].y << "," << vp->frc[deg].z << "\n";
    }
  }

  {
    std::ofstream ofs(f_lt);
    ofs << "step,line_id,lt\n";
    ofs << std::setprecision(17);
    for (int i = 0; i < (int)p_g->p_l.size(); i++) {
      const Line *lp = p_g->p_l[i];
      ofs << p_g->step << "," << i << "," << lp->lt << "\n";
    }
  }
}

} // namespace force_debug
