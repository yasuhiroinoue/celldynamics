/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "../restructure.h"
#include "../_class.h"
#include "../vec.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <set>
#include <cassert>
#include <cmath>

namespace restructure {

void removeTriangleVoid(Global *p_g, int vidx1, int vidx2, int vidx3) {
  return;
}

void removeTriangleVoid(Global *p_g) {
  for(int lidx = 0; lidx < (int)p_g->p_l.size(); ++lidx) {
    if(p_g->p_l[lidx]->ci.size() > 1) {
      continue;
    }

    int vidx1 = p_g->p_l[lidx]->vi[0];
    int vidx2 = p_g->p_l[lidx]->vi[1];

    for(int lidx2: p_g->p_v[vidx2]->li) {
      if(lidx2 <= lidx || p_g->p_l[lidx2]->ci.size() > 1) {
        continue;
      }

      int vidx3 = p_g->p_l[lidx2]->vi[0];
      if(vidx3 == vidx2) {
        vidx3 = p_g->p_l[lidx2]->vi[1];
      }

      for(int lidx3: p_g->p_v[vidx3]->li) {
        if(lidx3 <= lidx2 || p_g->p_l[lidx3]->ci.size() > 1) {
          continue;
        }

        if((p_g->p_l[lidx3]->vi[0] == vidx3 && p_g->p_l[lidx3]->vi[1] == vidx1) || (p_g->p_l[lidx3]->vi[0] == vidx1 && p_g->p_l[lidx3]->vi[1] == vidx3)) {
          std::cout << "remove triangle void " << vidx1 << " " << vidx2 << " " << vidx3 << std::endl;
          removeTriangleVoid(p_g, vidx1, vidx2, vidx3);
        }
      }
    }
  }
}

void killCellula(Global *p_g, int cidx) {
  for(int vidx: p_g->p_c[cidx]->vi) {
    findAndErase(p_g->p_v[vidx]->ci, cidx);
  }
  for(int lidx: p_g->p_c[cidx]->li) {
    findAndErase(p_g->p_l[lidx]->ci, cidx);
  }
  p_g->p_c[cidx]->vi.clear();
  p_g->p_c[cidx]->li.clear();
}

void apoptosis(Global *p_g, double time_constant) {
  for(int cidx: p_g->shrinking_ci) {
    p_g->p_c[cidx]->K_AREA /= pow(time_constant, 4);
    p_g->p_c[cidx]->AREA_EQ *= pow(time_constant, 2);
    for(int lidx: p_g->p_c[cidx]->li) {
      p_g->p_l[lidx]->K1_LENGTH /= time_constant * time_constant;
      //p_g->p_l[lidx]->K1_PCP_LENGTH /= time_constant;
      p_g->p_l[lidx]->K2_LENGTH /= time_constant * time_constant;
      p_g->p_l[lidx]->LENGTH_EQ *= time_constant;
    }
  }
}

} // namespace restructure
