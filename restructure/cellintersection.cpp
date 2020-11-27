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

bool containvc(Global *p_g, int vidx, int cidx) {
  bool ret = false;
  int csize = p_g->p_c[cidx]->vi.size();
  for(int i = 0; i < csize; ++i) {
    int tmp_vidx = p_g->p_c[cidx]->vi[i];
    int tmp_next_vidx = p_g->p_c[cidx]->vi[(i + 1) % csize];
    _vec<double> p1 = p_g->p_v[tmp_vidx]->loc[0] - p_g->p_v[vidx]->loc[0];
    _vec<double> p2 = p_g->p_v[tmp_next_vidx]->loc[0] - p_g->p_v[vidx]->loc[0];
    if(p1.y > p2.y) std::swap(p1, p2);
    if(p1.y <= 0 && 0 < p2.y) {
      if((p1 % p2).z < 0) ret = !ret;
    }
  }
  return ret;
}

_vec<double> projection(const pLine &l, const _vec<double> &p) {
  double t = ((p - l.first) * (l.first - l.second)) / (l.first - l.second).norm();
  _vec<double> ret = l.first + t * (l.first - l.second);
  return ret;
}

double distanceLP(const pLine &l, const _vec<double> &p) {
  return (p - projection(l, p)).norm();
}

int closestLineIndex(Global *p_g, int vidx, int cidx) {
  double dist = -1;
  int ret = 0;
  for(int lidx: p_g->p_c[cidx]->li) {
    int epoint_idx1 = p_g->p_l[lidx]->vi[0];
    int epoint_idx2 = p_g->p_l[lidx]->vi[1];
    pLine l = std::make_pair(p_g->p_v[epoint_idx1]->loc[0], p_g->p_v[epoint_idx2]->loc[0]);
    double tmp_dist = distanceLP(l, p_g->p_v[vidx]->loc[0]);
    if(dist == -1 || tmp_dist < dist) {
      dist = tmp_dist;
      ret = lidx;
    }
  }
  return ret;
}

void cellIntersectionType220(Global *p_g, int vinter_idx, int cobs_idx, int lobs_idx) {
  int vobs1_idx = p_g->p_l[lobs_idx]->vi[0];
  int vobs2_idx = p_g->p_l[lobs_idx]->vi[1];

  std::vector<std::pair<int,int>> epoint_vinter_idx;
  for(int lidx: p_g->p_v[vinter_idx]->li) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lidx]->vi[i] != vinter_idx) {
        epoint_vinter_idx.push_back(std::make_pair(p_g->p_l[lidx]->vi[i], lidx));
      }
    }
  }

  _vec<double> vinter = p_g->p_v[vinter_idx]->loc[0];
  _vec<double> vobs1 = p_g->p_v[vobs1_idx]->loc[0];

  int linter1_idx = -1;
  int linter2_idx = -1;
  int vadj1_idx = -1;
  int vadj2_idx = -1;
  _vec<double> vadj1_tmp = p_g->p_v[epoint_vinter_idx[0].first]->loc[0];
  _vec<double> vadj2_tmp = p_g->p_v[epoint_vinter_idx[1].first]->loc[0];
  if(((vadj1_tmp - vinter) % (vobs1 - vinter)).z * ((vadj1_tmp - vinter) % (vadj2_tmp - vinter)).z < 0) {
    linter1_idx = epoint_vinter_idx[0].second;
    linter2_idx = epoint_vinter_idx[1].second;
    vadj1_idx = epoint_vinter_idx[0].first;
    vadj2_idx = epoint_vinter_idx[1].first;
  }else{
    linter1_idx = epoint_vinter_idx[1].second;
    linter2_idx = epoint_vinter_idx[0].second;
    vadj1_idx = epoint_vinter_idx[1].first;
    vadj2_idx = epoint_vinter_idx[0].first;
  }

  int cinter_idx = p_g->p_l[linter1_idx]->ci[0];

  _vec<double> vadj1 = p_g->p_v[vadj1_idx]->loc[0];
  _vec<double> vadj2 = p_g->p_v[vadj2_idx]->loc[0];

  pLine lobs = std::make_pair(p_g->p_v[vobs1_idx]->loc[0], p_g->p_v[vobs2_idx]->loc[0]);
  pLine linter1 = std::make_pair(vinter, vadj1);
  pLine linter2 = std::make_pair(vinter, vadj2);

  int v1_idx = p_g->p_v.size();
  int l1_idx = p_g->p_l.size();
  int l2_idx = p_g->p_l.size() + 1;

  // Vertex v1
  Vertex *v1 = new Vertex;
  v1->li.push_back(linter2_idx);
  v1->li.push_back(l2_idx);
  v1->li.push_back(lobs_idx);
  v1->ci.push_back(cobs_idx);
  v1->ci.push_back(cinter_idx);
  v1->loc[0] = crossPoint(lobs, linter2);

  // Vertex vinter
  findAndErase(p_g->p_v[vinter_idx]->li, linter2_idx);
  p_g->p_v[vinter_idx]->li.push_back(l1_idx);
  p_g->p_v[vinter_idx]->li.push_back(l2_idx);
  p_g->p_v[vinter_idx]->ci.push_back(cobs_idx);
  p_g->p_v[vinter_idx]->loc[0] = crossPoint(lobs, linter1);

  // Vertex vobs1
  findAndErase(p_g->p_v[vobs1_idx]->li, lobs_idx);
  p_g->p_v[vobs1_idx]->li.push_back(l1_idx);

  // Vertex vobs2

  // Vertex vadj1

  // Vertex vadj2

  // Line l1
  Line *l1 = new Line;
  l1->vi[0] = vobs1_idx;
  l1->vi[1] = vinter_idx;
  l1->ci.push_back(cobs_idx);
  l1->lt = 0;
  l1->K1_LENGTH = K1_LENGTH;
  l1->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l1->K2_LENGTH = K2_LENGTH;
  l1->LENGTH_EQ = LENGTH_EQ;

  // Line l2
  Line *l2 = new Line;
  l2->vi[0] = vinter_idx;
  l2->vi[1] = v1_idx;
  l2->ci.push_back(cobs_idx);
  l2->ci.push_back(cinter_idx);
  l2->lt = 0;
  l2->K1_LENGTH = K1_LENGTH;
  l2->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l2->K2_LENGTH = K2_LENGTH;
  l2->LENGTH_EQ = LENGTH_EQ;


  // Line linter1

  // Line linter2
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[linter2_idx]->vi[i] == vinter_idx) {
      p_g->p_l[linter2_idx]->vi[i] = v1_idx;
    }
  }

  // Line lobs
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[lobs_idx]->vi[i] == vobs1_idx) {
      p_g->p_l[lobs_idx]->vi[i] = v1_idx;
    }
  }

  // Cellula cobs
  p_g->p_c[cobs_idx]->vi.push_back(vinter_idx);
  p_g->p_c[cobs_idx]->vi.push_back(v1_idx);
  p_g->p_c[cobs_idx]->li.push_back(l1_idx);
  p_g->p_c[cobs_idx]->li.push_back(l2_idx);

  // Cellula cinter
  p_g->p_c[cinter_idx]->vi.push_back(v1_idx);
  p_g->p_c[cinter_idx]->li.push_back(l2_idx);

  p_g->p_v.push_back(v1);
  p_g->p_l.push_back(l1);
  p_g->p_l.push_back(l2);

}

void cellIntersectionType321(Global *p_g, int vinter_idx, int cobs_idx, int lobs_idx) {
  std::vector<std::pair<int,int>> epoint_vinter_idx;
  for(int lidx: p_g->p_v[vinter_idx]->li) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lidx]->vi[i] != vinter_idx) {
        epoint_vinter_idx.push_back(std::make_pair(p_g->p_l[lidx]->vi[i], lidx));
      }
    }
  }

  int vobs1_idx = -1;
  int vobs2_idx = -1;
  int linter1_idx = -1;

  for(auto p: epoint_vinter_idx) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lobs_idx]->vi[i] == p.first) {
        vobs1_idx = p.first;
        vobs2_idx = p_g->p_l[lobs_idx]->vi[1 - i];
        linter1_idx = p.second;
      }
    }
  }
  assert(vobs1_idx != -1);
  assert(vobs2_idx != -1);
  assert(linter1_idx != -1);

  int linter2_idx = -1;
  int vadj1_idx = -1;
  for(auto p: epoint_vinter_idx) {
    if(p_g->p_l[p.second]->ci.size() == 2) {
      linter2_idx = p.second;
      vadj1_idx = p.first;
    }
  }
  assert(linter2_idx != -1);
  assert(vadj1_idx != -1);

  int linter3_idx = -1;
  int vadj2_idx = -1;
  for(auto p: epoint_vinter_idx) {
    if(p_g->p_l[p.second]->ci.size() == 1 && p.second != linter1_idx) {
      linter3_idx = p.second;
      vadj2_idx = p.first;
    }
  }
  assert(linter3_idx != -1);
  assert(vadj2_idx != -1);

  assert(p_g->p_l[linter1_idx]->ci.size() == 1);
  assert(p_g->p_l[linter3_idx]->ci.size() == 1);
  //int cinter1_idx = p_g->p_l[linter1_idx]->ci[0];
  int cinter2_idx = p_g->p_l[linter3_idx]->ci[0];

  int v1_idx = p_g->p_v.size();
  int l1_idx = p_g->p_l.size();

  pLine lobs = std::make_pair(p_g->p_v[vobs1_idx]->loc[0], p_g->p_v[vobs2_idx]->loc[0]);
  pLine linter2 = std::make_pair(p_g->p_v[vinter_idx]->loc[0], p_g->p_v[vadj1_idx]->loc[0]);
  pLine linter3 = std::make_pair(p_g->p_v[vinter_idx]->loc[0], p_g->p_v[vadj2_idx]->loc[0]);

  // Vertex v1
  Vertex *v1 = new Vertex;
  v1->li.push_back(lobs_idx);
  v1->li.push_back(l1_idx);
  v1->li.push_back(linter3_idx);
  v1->ci.push_back(cobs_idx);
  v1->ci.push_back(cinter2_idx);
  v1->loc[0] = crossPoint(lobs, linter3);

  // Vertex vobs1
  findAndErase(p_g->p_v[vobs1_idx]->li, lobs_idx);

  // Vertex vobs2

  // Vertex vinter
  findAndErase(p_g->p_v[vinter_idx]->li, linter3_idx);
  p_g->p_v[vinter_idx]->li.push_back(l1_idx);
  p_g->p_v[vinter_idx]->ci.push_back(cobs_idx);
  p_g->p_v[vinter_idx]->loc[0] = crossPoint(lobs, linter2);

  // Vertex vadj1

  // Vertex vadj2

  // Line l1
  Line *l1 = new Line;
  l1->vi[0] = v1_idx;
  l1->vi[1] = vinter_idx;
  l1->ci.push_back(cobs_idx);
  l1->ci.push_back(cinter2_idx);
  l1->lt = 0;
  l1->K1_LENGTH = K1_LENGTH;
  l1->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l1->K2_LENGTH = K2_LENGTH;
  l1->LENGTH_EQ = LENGTH_EQ;

  // Line lobs
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[lobs_idx]->vi[i] == vobs1_idx) {
      p_g->p_l[lobs_idx]->vi[i] = v1_idx;
    }
  }

  // Line linter1
  p_g->p_l[linter1_idx]->ci.push_back(cobs_idx);

  // Line linter2
  
  // Line linter3
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[linter3_idx]->vi[i] == vinter_idx) {
      p_g->p_l[linter3_idx]->vi[i] = v1_idx;
    }
  }

  // Cellula cobs
  p_g->p_c[cobs_idx]->vi.push_back(vinter_idx);
  p_g->p_c[cobs_idx]->vi.push_back(v1_idx);
  p_g->p_c[cobs_idx]->li.push_back(linter1_idx);
  p_g->p_c[cobs_idx]->li.push_back(l1_idx);

  // Cellula cinter1

  // Cellula cinter2
  p_g->p_c[cinter2_idx]->vi.push_back(v1_idx);
  p_g->p_c[cinter2_idx]->li.push_back(l1_idx);

  p_g->p_v.push_back(v1);
  p_g->p_l.push_back(l1);
}

void cellIntersectionType221(Global *p_g, int vinter_idx, int cobs_idx, int lobs_idx) {
  std::vector<std::pair<int,int>> epoint_vinter_idx;
  for(int lidx: p_g->p_v[vinter_idx]->li) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lidx]->vi[i] != vinter_idx) {
        epoint_vinter_idx.push_back(std::make_pair(p_g->p_l[lidx]->vi[i], lidx));
      }
    }
  }

  int vobs1_idx = -1;
  int vobs2_idx = -1;
  int linter1_idx = -1;

  for(auto p: epoint_vinter_idx) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lobs_idx]->vi[i] == p.first) {
        vobs1_idx = p.first;
        vobs2_idx = p_g->p_l[lobs_idx]->vi[1 - i];
        linter1_idx = p.second;
      }
    }
  }
  assert(vobs1_idx != -1);
  assert(vobs2_idx != -1);
  assert(linter1_idx != -1);

  int linter2_idx = -1;
  int vadj1_idx = -1;
  for(auto p: epoint_vinter_idx) {
    if(p.first != vobs1_idx) {
      linter2_idx = p.second;
      vadj1_idx = p.first;
    }
  }
  assert(linter2_idx != -1);
  assert(vadj1_idx != -1);

  assert(p_g->p_l[linter1_idx]->ci.size() == 1);
  //int cinter_idx = p_g->p_l[linter1_idx]->ci[0];

  pLine lobs = std::make_pair(p_g->p_v[vobs1_idx]->loc[0], p_g->p_v[vobs2_idx]->loc[0]);
  pLine linter2 = std::make_pair(p_g->p_v[vinter_idx]->loc[0], p_g->p_v[vadj1_idx]->loc[0]);

  // Vertex vinter
  p_g->p_v[vinter_idx]->li.push_back(lobs_idx);
  p_g->p_v[vinter_idx]->ci.push_back(cobs_idx);
  p_g->p_v[vinter_idx]->loc[0] = crossPoint(linter2, lobs);

  // Vertex vobs1
  findAndErase(p_g->p_v[vobs1_idx]->li, lobs_idx);

  // Vertex vobs2

  // Vertex vadj1

  // Line linter1
  p_g->p_l[linter1_idx]->ci.push_back(cobs_idx);

  // Line linter2

  // Line lobs
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[lobs_idx]->vi[i] == vobs1_idx) {
      p_g->p_l[lobs_idx]->vi[i] = vinter_idx;
    }
  }

  // Cellula cobs
  p_g->p_c[cobs_idx]->vi.push_back(vinter_idx);
  p_g->p_c[cobs_idx]->li.push_back(linter1_idx);

  // Cellula cinter
}

void cellIntersectionType330(Global *p_g, int vinter_idx, int cobs_idx, int lobs_idx) {
  int vobs1_idx = p_g->p_l[lobs_idx]->vi[0];
  int vobs2_idx = p_g->p_l[lobs_idx]->vi[1];

  std::vector<std::pair<int,int>> epoint_vinter_idx;
  for(int lidx: p_g->p_v[vinter_idx]->li) {
    for(int i = 0; i < 2; ++i) {
      if(p_g->p_l[lidx]->vi[i] != vinter_idx) {
        epoint_vinter_idx.push_back(std::make_pair(p_g->p_l[lidx]->vi[i], lidx));
      }
    }
  }

  int linter2_idx = -1;
  int vadj2_idx = -1;
  for(auto p: epoint_vinter_idx) {
    if(p_g->p_l[p.second]->ci.size() == 2) {
      linter2_idx = p.second;
      vadj2_idx = p.first;
    }
  }
  assert(linter2_idx != -1);
  assert(vadj2_idx != -1);

  _vec<double> vinter = p_g->p_v[vinter_idx]->loc[0];
  _vec<double> vobs1 = p_g->p_v[vobs1_idx]->loc[0];
  _vec<double> vadj2 = p_g->p_v[vadj2_idx]->loc[0];

  int linter1_idx = -1;
  int linter3_idx = -1;
  int vadj1_idx = -1;
  int vadj3_idx = -1;
  for(auto p: epoint_vinter_idx) {
    if(p_g->p_l[p.second]->ci.size() == 2) continue;
    _vec<double> vadj_tmp = p_g->p_v[p.first]->loc[0];
    if(((vadj2 - vinter) % (vadj_tmp - vinter)).z * ((vadj2 - vinter) % (vobs1 - vinter)).z > 0) {
      linter1_idx = p.second;
      vadj1_idx = p.first;
    }
    if(((vadj2 - vinter) % (vadj_tmp - vinter)).z * ((vadj2 - vinter) % (vobs1 - vinter)).z < 0) {
      linter3_idx = p.second;
      vadj3_idx = p.first;
    }
  }
  assert(linter1_idx != -1);
  assert(linter3_idx != -1);
  assert(vadj1_idx != -1);
  assert(vadj3_idx != -1);

  int cinter1_idx = p_g->p_l[linter1_idx]->ci[0];
  int cinter2_idx = p_g->p_l[linter3_idx]->ci[0];

  _vec<double> vadj1 = p_g->p_v[vadj1_idx]->loc[0];
  _vec<double> vadj3 = p_g->p_v[vadj3_idx]->loc[0];

  pLine lobs = std::make_pair(p_g->p_v[vobs1_idx]->loc[0], p_g->p_v[vobs2_idx]->loc[0]);
  pLine linter1 = std::make_pair(vinter, vadj1);
  pLine linter2 = std::make_pair(vinter, vadj2);
  pLine linter3 = std::make_pair(vinter, vadj3);

  int v1_idx = p_g->p_v.size();
  int v2_idx = p_g->p_v.size() + 1;
  int l1_idx = p_g->p_l.size();
  int l2_idx = p_g->p_l.size() + 1;
  int l3_idx = p_g->p_l.size() + 2;

  // Vertex v1
  Vertex *v1 = new Vertex;
  v1->li.push_back(linter2_idx);
  v1->li.push_back(l2_idx);
  v1->li.push_back(l3_idx);
  v1->ci.push_back(cobs_idx);
  v1->ci.push_back(cinter1_idx);
  v1->ci.push_back(cinter2_idx);
  v1->loc[0] = crossPoint(lobs, linter2);

  // Vertex v2
  Vertex *v2 = new Vertex;
  v2->li.push_back(linter3_idx);
  v2->li.push_back(l3_idx);
  v2->li.push_back(lobs_idx);
  v2->ci.push_back(cobs_idx);
  v2->ci.push_back(cinter2_idx);
  v2->loc[0] = crossPoint(lobs, linter3);

  // Vertex vinter
  findAndErase(p_g->p_v[vinter_idx]->li, linter2_idx);
  findAndErase(p_g->p_v[vinter_idx]->li, linter3_idx);
  p_g->p_v[vinter_idx]->li.push_back(l1_idx);
  p_g->p_v[vinter_idx]->li.push_back(l2_idx);
  findAndErase(p_g->p_v[vinter_idx]->ci, cinter2_idx);
  p_g->p_v[vinter_idx]->ci.push_back(cobs_idx);
  p_g->p_v[vinter_idx]->loc[0] = crossPoint(lobs, linter1);

  // Vertex vobs1
  findAndErase(p_g->p_v[vobs1_idx]->li, lobs_idx);
  p_g->p_v[vobs1_idx]->li.push_back(l1_idx);
  
  // Vertex vobs2

  // Vertex vadj1

  // Vertex vadj2

  // Vertex vadj3

  // Line l1
  Line *l1 = new Line;
  l1->vi[0] = vobs1_idx;
  l1->vi[1] = vinter_idx;
  l1->ci.push_back(cobs_idx);
  l1->lt = 0;
  l1->K1_LENGTH = K1_LENGTH;
  l1->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l1->K2_LENGTH = K2_LENGTH;
  l1->LENGTH_EQ = LENGTH_EQ;

  // Line l2
  Line *l2 = new Line;
  l2->vi[0] = vinter_idx;
  l2->vi[1] = v1_idx;
  l2->ci.push_back(cobs_idx);
  l2->ci.push_back(cinter1_idx);
  l2->lt = 0;
  l2->K1_LENGTH = K1_LENGTH;
  l2->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l2->K2_LENGTH = K2_LENGTH;
  l2->LENGTH_EQ = LENGTH_EQ;

  // Line l3
  Line *l3 = new Line;
  l3->vi[0] = v1_idx;
  l3->vi[1] = v2_idx;
  l3->ci.push_back(cobs_idx);
  l3->ci.push_back(cinter2_idx);
  l3->lt = 0;  
  l3->K1_LENGTH = K1_LENGTH;
  l3->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l3->K2_LENGTH = K2_LENGTH;
  l3->LENGTH_EQ = LENGTH_EQ;

  // Line linter1

  // Line linter2
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[linter2_idx]->vi[i] == vinter_idx) {
      p_g->p_l[linter2_idx]->vi[i] = v1_idx;
    }
  }

  // Line linter3
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[linter3_idx]->vi[i] == vinter_idx) {
      p_g->p_l[linter3_idx]->vi[i] = v2_idx;
    }
  }

  // Line lobs
  for(int i = 0; i < 2; ++i) {
    if(p_g->p_l[lobs_idx]->vi[i] == vobs1_idx) {
      p_g->p_l[lobs_idx]->vi[i] = v2_idx;
    }
  }

  // Cellula cobs
  p_g->p_c[cobs_idx]->vi.push_back(vinter_idx);
  p_g->p_c[cobs_idx]->vi.push_back(v1_idx);
  p_g->p_c[cobs_idx]->vi.push_back(v2_idx);
  p_g->p_c[cobs_idx]->li.push_back(l1_idx);
  p_g->p_c[cobs_idx]->li.push_back(l2_idx);
  p_g->p_c[cobs_idx]->li.push_back(l3_idx);

  // Cellula cinter1
  p_g->p_c[cinter1_idx]->vi.push_back(v1_idx);
  p_g->p_c[cinter1_idx]->li.push_back(l2_idx);

  // Cellula cinter2
  findAndErase(p_g->p_c[cinter2_idx]->vi, vinter_idx);
  p_g->p_c[cinter2_idx]->vi.push_back(v1_idx);
  p_g->p_c[cinter2_idx]->vi.push_back(v2_idx);
  p_g->p_c[cinter2_idx]->li.push_back(l3_idx);

  p_g->p_v.push_back(v1);
  p_g->p_v.push_back(v2);
  p_g->p_l.push_back(l1);
  p_g->p_l.push_back(l2);
  p_g->p_l.push_back(l3);
}

void cellIntersection(Global *p_g) {
  for(int vidx = 0; vidx < (int)p_g->p_v.size(); ++vidx) {
    for(int cidx = 0; cidx < (int)p_g->p_c.size(); ++cidx) {
      if(p_g->p_v[vidx]->ci.size() == 0) break;
      if(std::find(p_g->p_c[cidx]->vi.begin(), p_g->p_c[cidx]->vi.end(), vidx) != p_g->p_c[cidx]->vi.end()) {
        continue;
      }
      if(!containvc(p_g, vidx, cidx)) {
        continue;
      }

      int lobs_idx = closestLineIndex(p_g, vidx, cidx);

      if(p_g->p_l[lobs_idx]->ci.size() > 1) {
        continue;
      }

      int vinter_degree = p_g->p_v[vidx]->li.size();
      int vinter_outside = vinter_degree;
      int vinter_lobs = 0;

      std::vector<int> vadj_idx;
      for(int lidx: p_g->p_v[vidx]->li) {
        for(int i = 0; i < 2; ++i) {
          if(p_g->p_l[lidx]->vi[i] == vidx) continue;
          vadj_idx.push_back(p_g->p_l[lidx]->vi[i]);
        }
      }

      for(int i = 0; i < 2; ++i) {
        if(std::find(vadj_idx.begin(), vadj_idx.end(), p_g->p_l[lobs_idx]->vi[i]) != vadj_idx.end()) {
          --vinter_outside;
          ++vinter_lobs;
        }
      }
      for(int vv: p_g->p_c[cidx]->vi) {
        if(vv == p_g->p_l[lobs_idx]->vi[0] || vv == p_g->p_l[lobs_idx]->vi[1]) {
          continue;
        }
        if(std::find(vadj_idx.begin(), vadj_idx.end(), vv) != vadj_idx.end()) {
          --vinter_outside;
        }
      }

      std::cout << "  Intersection type: ";
      std::cout << vinter_degree << " " << vinter_outside << " " << vinter_lobs;
      std::cout << ", vidx: " << vidx << ", cidx: " << cidx << std::endl;

      bool dealed = true;
      if(vinter_degree == 2 && vinter_outside == 2 && vinter_lobs == 0) {
        cellIntersectionType220(p_g, vidx, cidx, lobs_idx);
        sortCounterClockwise(p_g);
      }else if(vinter_degree == 3 && vinter_outside == 2 && vinter_lobs == 1) {
        cellIntersectionType321(p_g, vidx, cidx, lobs_idx);
        sortCounterClockwise(p_g);
      }else if(vinter_degree == 2 && vinter_outside == 1 && vinter_lobs == 1) {
        cellIntersectionType221(p_g, vidx, cidx, lobs_idx);
        sortCounterClockwise(p_g);
      }else if(vinter_degree == 3 && vinter_outside == 3 && vinter_lobs == 0) {
        cellIntersectionType330(p_g, vidx, cidx, lobs_idx);
        sortCounterClockwise(p_g);
      }else{
        dealed = false;
      }

      if(!dealed) {
        std::cout << "unhandled" << std::endl;
      }else{
        std::cout << (isConsistent(p_g) ? "  consistent" : "anconsistent") << std::endl;
      }
    }
  }
}

}// namespace restructure
