/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "../restructure.h"
#include "../_class.h"
#include "../vec.h"
#include "../_parameters.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <cassert>
#include <cmath>

namespace restructure {

void cellRearrangeType442(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  int ladj1_idx = -1;
  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(lidx != lsep_idx && ladj1_idx == -1) {
      ladj1_idx = lidx;
    }else if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }

  int cadj1_idx = -1;
  int cadj2_idx = -1;
  int cadj4_idx = -1;
  for(int i = 0; i < 2; ++i) {
    int cidx = p_g->p_l[ladj1_idx]->ci[i];
    for(int j = 0; j < 2; ++j) {
      if(cidx != p_g->p_l[ladj2_idx]->ci[j]) {
        continue;
      }
      cadj1_idx = cidx;
      cadj2_idx = p_g->p_l[ladj2_idx]->ci[1 - j];
      cadj4_idx = p_g->p_l[ladj1_idx]->ci[1 - i];
    }
  }

  int cadj3_idx = -1;
  for(int cidx: p_g->p_v[vsep2_idx]->ci) {
    if(cidx != cadj2_idx && cidx != cadj4_idx) {
      cadj3_idx = cidx;
    }
  }

  int ladj3_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx == lsep_idx) {
      continue;
    }
    if(p_g->p_l[lidx]->ci[0] == cadj2_idx || p_g->p_l[lidx]->ci[1] == cadj2_idx) {
      ladj3_idx = lidx;
    }
  }

  int ladj4_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx && lidx != ladj3_idx) {
      ladj4_idx = lidx;
    }
  }

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, ladj1_idx);
  p_g->p_v[vsep1_idx]->li.push_back(ladj3_idx);
  findAndErase(p_g->p_v[vsep1_idx]->ci, cadj4_idx);
  p_g->p_v[vsep1_idx]->ci.push_back(cadj3_idx);

  p_g->p_v[vsep1_idx]->loc[0] = center;
  p_g->p_v[vsep1_idx]->loc[0].x -= (vsep1 - center).y;
  p_g->p_v[vsep1_idx]->loc[0].y += (vsep1 - center).x;

  // Vertex vsep2
  findAndErase(p_g->p_v[vsep2_idx]->li, ladj3_idx);
  p_g->p_v[vsep2_idx]->li.push_back(ladj1_idx);
  findAndErase(p_g->p_v[vsep2_idx]->ci, cadj2_idx);
  p_g->p_v[vsep2_idx]->ci.push_back(cadj1_idx);

  p_g->p_v[vsep2_idx]->loc[0] = center;
  p_g->p_v[vsep2_idx]->loc[0].x -= (vsep2 - center).y;
  p_g->p_v[vsep2_idx]->loc[0].y += (vsep2 - center).x;

  // Line lsep
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj2_idx);
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj4_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj1_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj3_idx);

  // Line ladj1
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj1_idx]->vi[i] == vsep1_idx) {
      p_g->p_l[ladj1_idx]->vi[i] = vsep2_idx;
    }
  }

  // Line ladj2

  // Line ladj3
   for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj3_idx]->vi[i] == vsep2_idx) {
      p_g->p_l[ladj3_idx]->vi[i] = vsep1_idx;
    }
  }

  // Line ladj4

  // Cellula cadj1
  p_g->p_c[cadj1_idx]->vi.push_back(vsep2_idx);
  p_g->p_c[cadj1_idx]->li.push_back(lsep_idx);

  // Cellula cadj2
  findAndErase(p_g->p_c[cadj2_idx]->vi, vsep2_idx);
  findAndErase(p_g->p_c[cadj2_idx]->li, lsep_idx);

  // Cellula cadj3
  p_g->p_c[cadj3_idx]->vi.push_back(vsep1_idx);
  p_g->p_c[cadj3_idx]->li.push_back(lsep_idx);

  // Cellula cadj4
  findAndErase(p_g->p_c[cadj4_idx]->vi, vsep1_idx);
  findAndErase(p_g->p_c[cadj4_idx]->li, lsep_idx);
}

void cellRearrangeType432(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  int ladj1_idx = -1;
  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(lidx != lsep_idx && ladj1_idx == -1) {
      ladj1_idx = lidx;
    }else if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }

  int cadj1_idx = -1;
  int cadj2_idx = -1;
  int cadj4_idx = -1;
  for(int i = 0; i < 2; ++i) {
    int cidx = p_g->p_l[ladj1_idx]->ci[i];
    for(int j = 0; j < 2; ++j) {
      if(cidx != p_g->p_l[ladj2_idx]->ci[j]) {
        continue;
      }
      cadj1_idx = cidx;
      cadj2_idx = p_g->p_l[ladj2_idx]->ci[1 - j];
      cadj4_idx = p_g->p_l[ladj1_idx]->ci[1 - i];
    }
  }

  int ladj3_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx && p_g->p_l[lidx]->ci[0] == cadj2_idx) {
      ladj3_idx = lidx;
    }
  }

  int ladj4_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx && lidx != ladj3_idx) {
      ladj4_idx = lidx;
    }
  }

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, ladj1_idx);
  p_g->p_v[vsep1_idx]->li.push_back(ladj3_idx);
  findAndErase(p_g->p_v[vsep1_idx]->ci, cadj4_idx);

  p_g->p_v[vsep1_idx]->loc[0] = center;
  p_g->p_v[vsep1_idx]->loc[0].x -= (vsep1 - center).y;
  p_g->p_v[vsep1_idx]->loc[0].y += (vsep1 - center).x;

  // Vertex vsep2
  findAndErase(p_g->p_v[vsep2_idx]->li, ladj3_idx);
  p_g->p_v[vsep2_idx]->li.push_back(ladj1_idx);
  findAndErase(p_g->p_v[vsep2_idx]->ci, cadj2_idx);
  p_g->p_v[vsep2_idx]->ci.push_back(cadj1_idx);

  p_g->p_v[vsep2_idx]->loc[0] = center;
  p_g->p_v[vsep2_idx]->loc[0].x -= (vsep2 - center).y;
  p_g->p_v[vsep2_idx]->loc[0].y += (vsep2 - center).x;

  // Line lsep
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj2_idx);
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj4_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj1_idx);

  // Line ladj1
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj1_idx]->vi[i] == vsep1_idx) {
      p_g->p_l[ladj1_idx]->vi[i] = vsep2_idx;
    }
  }

  // Line ladj2

  // Line ladj3
   for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj3_idx]->vi[i] == vsep2_idx) {
      p_g->p_l[ladj3_idx]->vi[i] = vsep1_idx;
    }
  }

  // Line ladj4

  // Cellula cadj1
  p_g->p_c[cadj1_idx]->vi.push_back(vsep2_idx);
  p_g->p_c[cadj1_idx]->li.push_back(lsep_idx);

  // Cellula cadj2
  findAndErase(p_g->p_c[cadj2_idx]->vi, vsep2_idx);
  findAndErase(p_g->p_c[cadj2_idx]->li, lsep_idx);

  // Cellula cadj4
  findAndErase(p_g->p_c[cadj4_idx]->vi, vsep1_idx);
  findAndErase(p_g->p_c[cadj4_idx]->li, lsep_idx);
}

void cellRearrangeType431(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  int ladj1_idx = -1;
  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(lidx != lsep_idx && p_g->p_l[lidx]->ci.size() == 2) {
      ladj1_idx = lidx;
    }else if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }
  assert(ladj1_idx != -1);
  assert(ladj2_idx != -1);

  int cadj1_idx = p_g->p_l[ladj2_idx]->ci[0];
  int cadj4_idx = p_g->p_l[lsep_idx]->ci[0];

  int cadj3_idx = -1;
  for(int cidx: p_g->p_v[vsep2_idx]->ci) {
    if(cidx != cadj4_idx) {
      cadj3_idx = cidx;
    }
  }

  int ladj3_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx && p_g->p_l[lidx]->ci.size() == 1) {
      ladj3_idx = lidx;
    }
  }

  int ladj4_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx && lidx != ladj3_idx) {
      ladj4_idx = lidx;
    }
  }

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, ladj1_idx);
  p_g->p_v[vsep1_idx]->li.push_back(ladj3_idx);
  findAndErase(p_g->p_v[vsep1_idx]->ci, cadj4_idx);
  p_g->p_v[vsep1_idx]->ci.push_back(cadj3_idx);

  p_g->p_v[vsep1_idx]->loc[0] = center;
  p_g->p_v[vsep1_idx]->loc[0].x -= (vsep1 - center).y;
  p_g->p_v[vsep1_idx]->loc[0].y += (vsep1 - center).x;

  // Vertex vsep2
  findAndErase(p_g->p_v[vsep2_idx]->li, ladj3_idx);
  p_g->p_v[vsep2_idx]->li.push_back(ladj1_idx);
  p_g->p_v[vsep2_idx]->ci.push_back(cadj1_idx);

  p_g->p_v[vsep2_idx]->loc[0] = center;
  p_g->p_v[vsep2_idx]->loc[0].x -= (vsep2 - center).y;
  p_g->p_v[vsep2_idx]->loc[0].y += (vsep2 - center).x;

  // Line lsep
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj4_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj1_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj3_idx);

  // Line ladj1
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj1_idx]->vi[i] == vsep1_idx) {
      p_g->p_l[ladj1_idx]->vi[i] = vsep2_idx;
    }
  }

  // Line ladj2

  // Line ladj3
   for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj3_idx]->vi[i] == vsep2_idx) {
      p_g->p_l[ladj3_idx]->vi[i] = vsep1_idx;
    }
  }

  // Line ladj4

  // Cellula cadj1
  p_g->p_c[cadj1_idx]->vi.push_back(vsep2_idx);
  p_g->p_c[cadj1_idx]->li.push_back(lsep_idx);

  // Cellula cadj3
  p_g->p_c[cadj3_idx]->vi.push_back(vsep1_idx);
  p_g->p_c[cadj3_idx]->li.push_back(lsep_idx);

  // Cellula cadj4
  findAndErase(p_g->p_c[cadj4_idx]->vi, vsep1_idx);
  findAndErase(p_g->p_c[cadj4_idx]->li, lsep_idx);
}

void cellRearrangeType332(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  int ladj1_idx = -1;
  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(lidx != lsep_idx && ladj1_idx == -1) {
      ladj1_idx = lidx;
    }else if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }

  int cadj1_idx = -1;
  int cadj2_idx = -1;
  int cadj3_idx = -1;
  for(int i = 0; i < 2; ++i) {
    int cidx = p_g->p_l[ladj1_idx]->ci[i];
    for(int j = 0; j < 2; ++j) {
      if(cidx != p_g->p_l[ladj2_idx]->ci[j]) {
        continue;
      }
      cadj1_idx = cidx;
      cadj2_idx = p_g->p_l[ladj2_idx]->ci[1 - j];
      cadj3_idx = p_g->p_l[ladj1_idx]->ci[1 - i];
    }
  }

  int ladj3_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx) {
      ladj3_idx = lidx;
    }
  }

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, ladj1_idx);
  p_g->p_v[vsep1_idx]->li.push_back(ladj3_idx);

  p_g->p_v[vsep1_idx]->loc[0] = center;
  p_g->p_v[vsep1_idx]->loc[0].x -= (vsep1 - center).y;
  p_g->p_v[vsep1_idx]->loc[0].y += (vsep1 - center).x;

  // Vertex vsep2
  findAndErase(p_g->p_v[vsep2_idx]->li, ladj3_idx);
  p_g->p_v[vsep2_idx]->li.push_back(ladj1_idx);
  findAndErase(p_g->p_v[vsep2_idx]->ci, cadj2_idx);
  p_g->p_v[vsep2_idx]->ci.push_back(cadj1_idx);

  p_g->p_v[vsep2_idx]->loc[0] = center;
  p_g->p_v[vsep2_idx]->loc[0].x -= (vsep2 - center).y;
  p_g->p_v[vsep2_idx]->loc[0].y += (vsep2 - center).x;

  // Line lsep
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj2_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj1_idx);

  // Line ladj1
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj1_idx]->vi[i] == vsep1_idx) {
      p_g->p_l[ladj1_idx]->vi[i] = vsep2_idx;
    }
  }

  // Line ladj2

  // Line ladj3
   for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj3_idx]->vi[i] == vsep2_idx) {
      p_g->p_l[ladj3_idx]->vi[i] = vsep1_idx;
    }
  }

  // Cellula cadj1
  p_g->p_c[cadj1_idx]->vi.push_back(vsep2_idx);
  p_g->p_c[cadj1_idx]->li.push_back(lsep_idx);

  // Cellula cadj2
  findAndErase(p_g->p_c[cadj2_idx]->vi, vsep2_idx);
  findAndErase(p_g->p_c[cadj2_idx]->li, lsep_idx);

  // Cellula cadj3
}

void cellRearrangeType321(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  int ladj1_idx = -1;
  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(p_g->p_l[lidx]->ci.size() == 2) {
      ladj1_idx = lidx;
    }else if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }

  int cadj1_idx = p_g->p_l[ladj2_idx]->ci[0];
  int cadj4_idx = p_g->p_l[lsep_idx]->ci[0];

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, ladj1_idx);
  findAndErase(p_g->p_v[vsep1_idx]->ci, cadj4_idx);

  p_g->p_v[vsep1_idx]->loc[0] = center;
  p_g->p_v[vsep1_idx]->loc[0].x -= (vsep1 - center).y;
  p_g->p_v[vsep1_idx]->loc[0].y += (vsep1 - center).x;

  // Vertex vsep2
  p_g->p_v[vsep2_idx]->li.push_back(ladj1_idx);
  p_g->p_v[vsep2_idx]->ci.push_back(cadj1_idx);

  p_g->p_v[vsep2_idx]->loc[0] = center;
  p_g->p_v[vsep2_idx]->loc[0].x -= (vsep2 - center).y;
  p_g->p_v[vsep2_idx]->loc[0].y += (vsep2 - center).x;

  // Line lsep
  findAndErase(p_g->p_l[lsep_idx]->ci, cadj4_idx);
  p_g->p_l[lsep_idx]->ci.push_back(cadj1_idx);

  // Line ladj1
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj1_idx]->vi[i] == vsep1_idx) {
      p_g->p_l[ladj1_idx]->vi[i] = vsep2_idx;
    }
  }

  // Line ladj2

  // Cellula cadj1
  p_g->p_c[cadj1_idx]->vi.push_back(vsep2_idx);
  p_g->p_c[cadj1_idx]->li.push_back(lsep_idx);

  // Cellula cadj4
  findAndErase(p_g->p_c[cadj4_idx]->vi, vsep1_idx);
  findAndErase(p_g->p_c[cadj4_idx]->li, lsep_idx);
}

void cellRearrangeType211(Global *p_g, int lsep_idx) {
  int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
  int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

  /*
  int ladj1_idx = -1;
  for(int lidx: p_g->p_v[vsep1_idx]->li) {
    if(lidx != lsep_idx) {
      ladj1_idx = lidx;
    }
  }
  */

  int ladj2_idx = -1;
  for(int lidx: p_g->p_v[vsep2_idx]->li) {
    if(lidx != lsep_idx) {
      ladj2_idx = lidx;
    }
  }

  int cadj = p_g->p_l[lsep_idx]->ci[0];

  _vec<double> vsep1 = p_g->p_v[vsep1_idx]->loc[0];
  _vec<double> vsep2 = p_g->p_v[vsep2_idx]->loc[0];
  _vec<double> center = 0.5 * (vsep1 + vsep2);

  // Vertex vsep1
  findAndErase(p_g->p_v[vsep1_idx]->li, lsep_idx);
  p_g->p_v[vsep1_idx]->li.push_back(ladj2_idx);
  p_g->p_v[vsep1_idx]->loc[0] = center;

  // Vertex vsep2
  p_g->p_v[vsep2_idx]->li.clear();
  p_g->p_v[vsep2_idx]->ci.clear();

  // Line lsep
  p_g->p_l[lsep_idx]->ci.clear();
  p_g->p_l[lsep_idx]->vi[1] = p_g->p_l[lsep_idx]->vi[0];

  // Line ladj1

  // Line ladj2
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[ladj2_idx]->vi[i] == vsep2_idx) {
      p_g->p_l[ladj2_idx]->vi[i] = vsep1_idx;
    }
  }

  // Cellula cadj
  findAndErase(p_g->p_c[cadj]->vi, vsep2_idx);
  findAndErase(p_g->p_c[cadj]->li, lsep_idx);
}

void cellRearrange2(Global *p_g) {
  for(int lsep_idx = 0; lsep_idx < (int)p_g->p_l.size(); ++lsep_idx) {
    int vsep1_idx = p_g->p_l[lsep_idx]->vi[0];
    int vsep2_idx = p_g->p_l[lsep_idx]->vi[1];

    if(p_g->p_l[lsep_idx]->ci.size() == 0) {
      continue;
    }

    if((p_g->p_v[vsep1_idx]->loc[0] - p_g->p_v[vsep2_idx]->loc[0]).norm() > K1_PCP_LENGTH) {
      continue;
    }

    bool triangle = false;
    for(int cidx: p_g->p_l[lsep_idx]->ci) {
      if(p_g->p_c[cidx]->vi.size() <= 3) {
        triangle = true;
      }
    }
    if(triangle) {
      continue;
    }

    if(p_g->p_v[vsep1_idx]->li.size() < p_g->p_v[vsep2_idx]->li.size()) {
      std::swap(vsep1_idx, vsep2_idx);
      std::swap(p_g->p_l[lsep_idx]->vi[0], p_g->p_l[lsep_idx]->vi[1]);
    }else if(p_g->p_v[vsep1_idx]->li.size() == p_g->p_v[vsep2_idx]->li.size()) {
      if(p_g->p_v[vsep1_idx]->ci.size() < p_g->p_v[vsep2_idx]->ci.size()) {
        std::swap(vsep1_idx, vsep2_idx);
        std::swap(p_g->p_l[lsep_idx]->vi[0], p_g->p_l[lsep_idx]->vi[1]);
      }
    }

    int line_cnt = p_g->p_v[vsep1_idx]->li.size() + p_g->p_v[vsep2_idx]->li.size() - 2;
    int cellula_total_cnt = p_g->p_v[vsep1_idx]->ci.size() + p_g->p_v[vsep2_idx]->ci.size() - p_g->p_l[lsep_idx]->ci.size();
    int cellula_lsep_cnt = p_g->p_l[lsep_idx]->ci.size();

    if(line_cnt == 4 && cellula_total_cnt == 4 && cellula_lsep_cnt == 2) {
      std::cout << "  rearrange(type442): Line-" << lsep_idx << std::endl;
      cellRearrangeType442(p_g, lsep_idx);
      //std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }
    else if(line_cnt == 4 && cellula_total_cnt == 3 && cellula_lsep_cnt == 2) {
      std::cout << "  rearrange(type432): Line-" << lsep_idx << std::endl;
      cellRearrangeType432(p_g, lsep_idx);
      //std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }else if(line_cnt == 4 && cellula_total_cnt == 3 && cellula_lsep_cnt == 1) {
      std::cout << "  rearrange(type431): Line-" << lsep_idx << std::endl;
      cellRearrangeType431(p_g, lsep_idx);
      //std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }else if(line_cnt == 3 && cellula_total_cnt == 3 && cellula_lsep_cnt == 2) {
      std::cout << "  rearrange(type332): Line-" << lsep_idx << std::endl;
      cellRearrangeType332(p_g, lsep_idx);
      //std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }else if(line_cnt == 3 && cellula_total_cnt == 2 && cellula_lsep_cnt == 1) {
      std::cout << "  rearrange(type321): Line-" << lsep_idx << std::endl;
      cellRearrangeType321(p_g, lsep_idx);
      //std::cout << (isConsistent(p_g) ? "consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }else if(line_cnt == 2 && cellula_total_cnt == 1 && cellula_lsep_cnt == 1){
      std::cout << "  rearrange(type211): Line-" << lsep_idx << std::endl;
      cellRearrangeType211(p_g, lsep_idx);
      std::cout << (isConsistent(p_g) ? "  consistent" : "anconsistent") << std::endl;
      sortCounterClockwise(p_g);
    }else{
      std::cout << "unhandled rearrange(type" << line_cnt << cellula_total_cnt << cellula_lsep_cnt << ") Line-" << lsep_idx << std::endl;
      std::cout << "vidx: " << p_g->p_l[lsep_idx]->vi[0] << " " << p_g->p_l[lsep_idx]->vi[1] << std::endl;
    }
  }
}

} //namespace restructure
