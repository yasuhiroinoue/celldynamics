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
#include <cassert>
#include <cmath>

namespace restructure {

void vertexOverlap(Global *p_g) {
  for(int lol_idx = 0; lol_idx < (int)p_g->p_l.size(); ++lol_idx) {
    if(p_g->p_l[lol_idx]->ci.size() == 0) continue;
    int vol1_idx = p_g->p_l[lol_idx]->vi[0];
    int vol2_idx = p_g->p_l[lol_idx]->vi[1];
    
    if((p_g->p_v[vol1_idx]->loc[0] - p_g->p_v[vol2_idx]->loc[0]).norm() > 1e-5) continue;

    std::cout << "start to unite vertexs" << std::endl;
    for(int lidx: p_g->p_v[vol2_idx]->li) {
      for(int i = 0; i < 2; ++i) {
        if(p_g->p_l[lidx]->vi[i] == vol2_idx) {
          p_g->p_l[lidx]->vi[i] = vol1_idx;
        }
      }
    }
    for(int cidx: p_g->p_v[vol2_idx]->ci) {
      for(int i = 0; i < (int)p_g->p_c[cidx]->vi.size(); ++i) {
        if(p_g->p_c[cidx]->vi[i] == vol2_idx) {
          findAndErase(p_g->p_c[cidx]->vi, vol1_idx);
          p_g->p_c[cidx]->vi[i] = vol1_idx;
        }
      }
    }

    for(int lidx: p_g->p_v[vol2_idx]->li) {
      if(std::find(p_g->p_v[vol1_idx]->li.begin(), p_g->p_v[vol1_idx]->li.end(), lidx) == p_g->p_v[vol1_idx]->li.end()) {
        p_g->p_v[vol1_idx]->li.push_back(lidx);
      }
    }
    for(int cidx: p_g->p_v[vol2_idx]->ci) {
      if(std::find(p_g->p_v[vol1_idx]->ci.begin(), p_g->p_v[vol1_idx]->ci.end(), cidx) == p_g->p_v[vol1_idx]->ci.end()) {
        p_g->p_v[vol1_idx]->ci.push_back(cidx);
      }
    }

    for(int cidx: p_g->p_l[lol_idx]->ci) {
      findAndErase(p_g->p_c[cidx]->li, lol_idx);
    }

    p_g->p_v[vol2_idx]->li.clear();
    p_g->p_v[vol2_idx]->ci.clear();
    p_g->p_v[vol2_idx]->loc[0] = _vec<double>(0, 0, 0);
    p_g->p_l[lol_idx]->ci.clear();
    for(int i = 0; i < 2; ++i) {
      p_g->p_l[lol_idx]->vi[i] = vol2_idx;
    }

    std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
    sortCounterClockwise(p_g);
    std::cout << "Unite Vertex-" << vol1_idx << " and Vertex-" << vol2_idx << std::endl;
  }
}

} // namespace restructure
