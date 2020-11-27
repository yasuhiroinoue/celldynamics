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

bool isConsistent(const Global *p_g) {
  for(int vidx = 0; vidx < (int)p_g->p_v.size(); ++vidx) {
    if(p_g->p_v[vidx]->ci.size() == 0) continue;
    for(int lidx: p_g->p_v[vidx]->li) {
      std::vector<int> ll = {p_g->p_l[lidx]->vi[0], p_g->p_l[lidx]->vi[1]};
      if(find(ll.begin(), ll.end(), vidx) == ll.end()) {
        std::cout << "Line-" << lidx << " should have vertex-" << vidx << std::endl;
        return false;
      }
    }
    for(int cidx: p_g->p_v[vidx]->ci) {
      std::vector<int> cc = p_g->p_c[cidx]->vi;
      if(find(cc.begin(), cc.end(), vidx) == cc.end()) {
        std::cout << "Cellula-" << cidx << "should have vertex-" << vidx << std::endl;
      }
    }
  }

  for(int lidx = 0; lidx < (int)p_g->p_l.size(); ++lidx) {
    if(p_g->p_l[lidx]->ci.size() == 0) continue;
    for(int i = 0; i < 2; ++i) {
      std::vector<int> vv = p_g->p_v[p_g->p_l[lidx]->vi[i]]->li;
      if(find(vv.begin(), vv.end(), lidx) == vv.end()) {
        std::cout << "Vertex-" << p_g->p_l[lidx]->vi[i] << "should have line-" << lidx << std::endl;
      }
    }
    for(int cidx: p_g->p_l[lidx]->ci) {
      std::vector<int> cc = p_g->p_c[cidx]->li;
      if(find(cc.begin(), cc.end(), lidx) == cc.end()) {
        std::cout << "Cellula-" << cidx << " should have line-" << lidx << std::endl;
      }
    }
  }

  for(int cidx = 0; cidx < (int)p_g->p_c.size(); ++cidx) {
    for(int vidx: p_g->p_c[cidx]->vi) {
      std::vector<int> vv = p_g->p_v[vidx]->ci;
      if(find(vv.begin(), vv.end(), cidx) == vv.end()) {
        std::cout << "Vertex-" << vidx << " should have cellula-" << cidx << std::endl;
      }
    }
    for(int lidx: p_g->p_c[cidx]->li) {
      std::vector<int> ll = p_g->p_l[lidx]->ci;
      if(find(ll.begin(), ll.end(), cidx) == ll.end()) {
        std::cout << "Line-" << lidx << " should have cellula-" << cidx << std::endl;
      }
    }
  }
  return true;
}

void sortCounterClockwise(Global *p_g) {
  // For each cellula, sort Cellula::vi and Cellula::li such that Cellula::vi is 
  // counterclockwise and Cellula::li[i] consists of Cellula::vi[i] and Cellula::vi[i + 1].
  // It is not guaranteed that p_l[Cellula::li[i]]->vi[0] is Cellula::vi[i].
  for(int cidx = 0; cidx < (int)p_g->p_c.size(); ++cidx) {
    std::vector<int> tmp_vi;
    std::vector<int> tmp_li;

    int prev_vidx = -1, curr_vidx = p_g->p_c[cidx]->vi[0];
    tmp_vi.push_back(curr_vidx);
    while(tmp_li.size() < p_g->p_c[cidx]->li.size()) {
      for(int lidx: p_g->p_c[cidx]->li) {
        if(p_g->p_l[lidx]->vi[0] == curr_vidx && p_g->p_l[lidx]->vi[1] != prev_vidx) {
          prev_vidx = curr_vidx;
          curr_vidx = p_g->p_l[lidx]->vi[1];
          tmp_li.push_back(lidx);
          if(curr_vidx != tmp_vi[0]){
            tmp_vi.push_back(curr_vidx);
          }
        }else if(p_g->p_l[lidx]->vi[0] != prev_vidx && p_g->p_l[lidx]->vi[1] == curr_vidx) {
          prev_vidx = curr_vidx;
          curr_vidx = p_g->p_l[lidx]->vi[0];
          tmp_li.push_back(lidx);
          if(curr_vidx != tmp_vi[0]){
            tmp_vi.push_back(curr_vidx);
          }
        }
      }
    }

    assert(tmp_vi.size() == p_g->p_c[cidx]->vi.size());
    assert(tmp_li.size() == p_g->p_c[cidx]->li.size());

    double area = 0.0;
    for(int i = 0; i < (int)tmp_vi.size(); ++i) {
      _vec<double> r1 = p_g->p_v[tmp_vi[i]]->loc[0];
      _vec<double> r2 = p_g->p_v[tmp_vi[(i + 1) % tmp_vi.size()]]->loc[0];
      area += 0.5 * (r1 % r2).z;
    }
    if(area < 0) {
      std::reverse(tmp_vi.begin(), tmp_vi.end());
      std::reverse(tmp_li.begin(), tmp_li.end());
    }

    for(int i = 0; i < (int)tmp_vi.size(); ++i) {
      p_g->p_c[cidx]->vi[i] = tmp_vi[i];
      p_g->p_c[cidx]->li[i] = tmp_li[i];
    }
  }
}

int ccw(_vec<double> a, _vec<double> _b, _vec<double> _c) {
  _vec<double> b = _b - a;
  _vec<double> c = _c - a;
  if((b % c).z > 0) return +1;  // counter clockwise
  if((b % c).z < 0) return -1;  // clockwise
  if(b * c < 0) return +2;  // c--a--b on line
  if(b.norm() < c.norm()) return -2;  // a--b--c on line
  return 0;
}

bool isConvex(const Global *p_g, int cidx) {
  int vsize = p_g->p_c[cidx]->vi.size();
  for(int i = 0; i < vsize; ++i) {
    _vec<double> prev = p_g->p_v[p_g->p_c[cidx]->vi[(i - 1 + vsize) % vsize]]->loc[0];
    _vec<double> curr = p_g->p_v[p_g->p_c[cidx]->vi[i]]->loc[0];
    _vec<double> next = p_g->p_v[p_g->p_c[cidx]->vi[(i + 1) % vsize]]->loc[0];
    if(ccw(prev, curr, next) < 0) return false;
  }
  return true;
}

double EPS = 1e-9;

bool isIntersectSegmentLine(pLine &s, pLine &l) {
  if(abs(((s.second - s.first) % (l.second - l.first)).z) < EPS) {
    return false;
  }
  _vec<double> m1 = (l.second - l.first) % (s.first - l.first);
  _vec<double> m2 = (l.second - l.first) % (s.second - l.first);
  if(m1.z * m2.z > EPS) return false;
  return true;
}

_vec<double> crossPoint(pLine &s, pLine &t) {
  _vec<double> sv = s.second - s.first;
  _vec<double> tv = t.second - t.first;
  assert((sv % tv).norm() > EPS);
  double length = (tv % (t.first - s.first)).z / (tv % sv).z;
  return s.first + sv * length;
}

void cellDivision(Global *p_g, int cellula_idx, _vec<double> axis) {
  Cellula *cp = p_g->p_c[cellula_idx];
  
  _vec<double> center = cp->center;
  pLine division = std::make_pair(center, center + axis);
  std::vector<std::pair<int, _vec<double>>> crosspoint;
  std::vector<int> epoint_idx;

  // sort lines counterclockwise
  // at the same time, sort p_g->p_l[lidx]->vi counterclockwise
  std::vector<int> anticlockwise_lidx;
  int cp_vsize = cp->vi.size();
  for(int i = 0; i < cp_vsize; ++i) {
    for(int lidx: cp->li) {
      if(p_g->p_l[lidx]->vi[0] == cp->vi[i] && p_g->p_l[lidx]->vi[1] == cp->vi[(i + 1) % cp_vsize]){
        anticlockwise_lidx.push_back(lidx);
        break;
      }else if(p_g->p_l[lidx]->vi[0] == cp->vi[(i + 1) % cp_vsize] && p_g->p_l[lidx]->vi[1] == cp->vi[i]) {
        std::swap(p_g->p_l[lidx]->vi[0], p_g->p_l[lidx]->vi[1]);
        anticlockwise_lidx.push_back(lidx);
      }
    }
  }

  for(int lidx: anticlockwise_lidx) {
    int tmp_epoint_idx1 = p_g->p_l[lidx]->vi[0];
    int tmp_epoint_idx2 = p_g->p_l[lidx]->vi[1];

    pLine segment = std::make_pair(*(p_g->p_v[tmp_epoint_idx1]->loc), *(p_g->p_v[tmp_epoint_idx2]->loc));
    if(isIntersectSegmentLine(segment, division)) {
      epoint_idx.push_back(tmp_epoint_idx1);
      epoint_idx.push_back(tmp_epoint_idx2);

      crosspoint.push_back(std::make_pair(lidx, crossPoint(segment, division)));
    }
  }

  if(crosspoint.size() != 2) {
    throw "Cannot divide cellulas";
    return;
  }

  int lim_num1 = -1, lim_num2 = -1;

  for(int i = 0; i < cp_vsize; ++i) {
    int vidx = cp->vi[i];
    if(vidx == epoint_idx[1]) {
      lim_num1 = i;
    }
    if(vidx == epoint_idx[2]) {
      lim_num2 = i - 1;
    }
  }

  assert(lim_num1 != -1);
  assert(lim_num2 != -1);

  int v1_idx = p_g->p_v.size();
  int v2_idx = p_g->p_v.size() + 1;

  int l1_idx = p_g->p_l.size();
  int l2_idx = p_g->p_l.size() + 1;
  int l3_idx = p_g->p_l.size() + 2;

  int c1_idx = p_g->p_c.size();

  // vertex v1
  Vertex *v1 = new Vertex;
  v1->li.push_back(crosspoint[0].first);
  v1->li.push_back(l1_idx);
  v1->li.push_back(l3_idx);
  for(int c: p_g->p_l[crosspoint[0].first]->ci) {
    v1->ci.push_back(c);
  }
  v1->ci.push_back(c1_idx);
  v1->loc[0] = crosspoint[0].second;
  v1->loc[0].z = 0.0;
  //v1->loc[1] = ???

  // vertex v2
  Vertex *v2 = new Vertex;
  v2->li.push_back(crosspoint[1].first);
  v2->li.push_back(l2_idx);
  v2->li.push_back(l3_idx);
  for(int c: p_g->p_l[crosspoint[1].first]->ci) {
    v2->ci.push_back(c);
  }
  v2->ci.push_back(c1_idx);
  v2->loc[0] = crosspoint[1].second;
  v2->loc[0].z = 0.0;
  //v2->loc[1] = ???

  // vertex epoint_idx[1]
  p_g->p_v[epoint_idx[1]]->li.push_back(l1_idx);
  findAndErase(p_g->p_v[epoint_idx[1]]->li, crosspoint[0].first);
  p_g->p_v[epoint_idx[1]]->ci.push_back(c1_idx);
  findAndErase(p_g->p_v[epoint_idx[1]]->ci, cellula_idx);

  // vertex epoint_idx[2]
  p_g->p_v[epoint_idx[2]]->li.push_back(l2_idx);
  findAndErase(p_g->p_v[epoint_idx[2]]->li, crosspoint[1].first);
  p_g->p_v[epoint_idx[2]]->ci.push_back(c1_idx);
  findAndErase(p_g->p_v[epoint_idx[2]]->ci, cellula_idx);

  for(int i = lim_num1 + 1; i <= lim_num2; ++i) {
    int vidx = cp->vi[i];
    p_g->p_v[vidx]->ci.push_back(c1_idx);
    findAndErase(p_g->p_v[vidx]->ci, cellula_idx);
  }

  // line l1
  Line *l1 = new Line;
  l1->vi[0] = v1_idx;
  l1->vi[1] = epoint_idx[1];
  for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
    l1->ci.push_back(cidx);
  }
  l1->ci.push_back(c1_idx);
  findAndErase(l1->ci, cellula_idx);
  l1->lt = 0.0;
  l1->K1_LENGTH = K1_LENGTH;
  l1->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l1->K2_LENGTH = K2_LENGTH;
  l1->LENGTH_EQ = LENGTH_EQ;

  // line l2
  Line *l2 = new Line;
  l2->vi[0] = v2_idx;
  l2->vi[1] = epoint_idx[2];
  for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
    l2->ci.push_back(cidx);
  }
  l2->ci.push_back(c1_idx);
  findAndErase(l2->ci, cellula_idx);
  l2->lt = 0.0;
  l2->K1_LENGTH = K1_LENGTH;
  l2->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l2->K2_LENGTH = K2_LENGTH;
  l2->LENGTH_EQ = LENGTH_EQ;

  // line l3
  Line *l3 = new Line;
  l3->vi[0] = v1_idx;
  l3->vi[1] = v2_idx;
  l3->ci.push_back(cellula_idx);
  l3->ci.push_back(c1_idx);
  l3->lt = 0.0;
  l3->K1_LENGTH = K1_LENGTH;
  l3->K1_PCP_LENGTH = K1_PCP_LENGTH;
  l3->K2_LENGTH = K2_LENGTH;
  l3->LENGTH_EQ = LENGTH_EQ;

  // line crosspoint[0].first
  p_g->p_l[crosspoint[0].first]->vi[1] = v1_idx;

  // line crosspoint[1].first
  p_g->p_l[crosspoint[1].first]->vi[0] = v2_idx;

  // line lim_num1-th --- lim_num2-th
  for(int i = lim_num1; i <= lim_num2; ++i) {
    int lidx = anticlockwise_lidx[i];
    p_g->p_l[lidx]->ci.push_back(c1_idx);
    findAndErase(p_g->p_l[lidx]->ci, cellula_idx);
  }

  // cellula c1
  Cellula *c1 = new Cellula;
  c1->vi.push_back(v1_idx);
  for(int i = lim_num1; i <= lim_num2 + 1; ++i) {
    int vidx = cp->vi[i];
    c1->vi.push_back(vidx);
  }
  c1->vi.push_back(v2_idx);
  c1->li.push_back(l1_idx);
  for(int i = lim_num1; i <= lim_num2; ++i) {
    int lidx = anticlockwise_lidx[i];
    c1->li.push_back(lidx);
  }
  c1->li.push_back(l2_idx);
  c1->li.push_back(l3_idx);
  c1->center = cp->center;
  c1->cell_time = 0.0;
  c1->cell_phase = 0.0; /////////////////////////
  c1->cell_T = Pulse_T;
  c1->fix = 0;
  c1->AREA_EQ = AREA_EQ;
  c1->K_AREA = K_AREA;

  // cellula cellula_idx
  cp->vi.erase(cp->vi.begin() + lim_num1, cp->vi.begin() + lim_num2 + 2);
  cp->vi.push_back(v1_idx);
  cp->vi.push_back(v2_idx);
  for(int i = lim_num1; i <= lim_num2; ++i) {
    //std::cerr << anticlockwise_lidx[i] << std::endl;
    findAndErase(cp->li, anticlockwise_lidx[i]);
  }
  cp->li.push_back(l3_idx);

  // cellula adjusting line-crosspoint[0].first (if exist)
  int adj_cidx1 = -1;
  for(int cidx: p_g->p_l[crosspoint[0].first]->ci) {
    if(cidx != cellula_idx) {
      adj_cidx1 = cidx;
    }
  }
  if(adj_cidx1 != -1) {
    p_g->p_c[adj_cidx1]->vi.push_back(v1_idx);
    p_g->p_c[adj_cidx1]->li.push_back(l1_idx);
  }

  // cellula adjusting line-crosspoint[1].first (if exist)
  int adj_cidx2 = -1;
  for(int cidx: p_g->p_l[crosspoint[1].first]->ci) {
    if(cidx != cellula_idx) {
      adj_cidx2 = cidx;
    }
  }
  if(adj_cidx2 != -1) {
    p_g->p_c[adj_cidx2]->vi.push_back(v2_idx);
    p_g->p_c[adj_cidx2]->li.push_back(l2_idx);
  }

  p_g->p_c.push_back(c1);
  p_g->p_l.push_back(l1);
  p_g->p_l.push_back(l2);
  p_g->p_l.push_back(l3);
  p_g->p_v.push_back(v1);
  p_g->p_v.push_back(v2);

  //isConsistent(p_g);

  //sortAntiClockwise(p_g);
  sortCounterClockwise(p_g);
}

}// namespace restructure
