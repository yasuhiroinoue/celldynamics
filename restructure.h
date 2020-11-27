/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _CRA_H
#define _CRA_H

#include "_class.h"
#include "vec.h"
#include <vector>
#include <iostream>
#include <algorithm>

namespace restructure {
  typedef std::pair<_vec<double>, _vec<double>> pLine;
  
  void sortAntiClockwise(Global *);
  void cellRearrange(Global *);
  void cellRearrange2(Global *);
  bool isConsistent(const Global *);
  bool isConvex(const Global * , int);
  _vec<double> crossPoint(pLine&, pLine&);
  void sortCounterClockwise(Global *);
  void cellDivision(Global *, int, _vec<double>);
  void cellIntersection(Global *);
  void vertexOverlap(Global *);
  void removeTriangleVoid(Global *);
  void apoptosis(Global *, double);

  template<typename T>
  void findAndErase(std::vector<T> &vector, T search) {
    auto itr = std::find(vector.begin(), vector.end(), search);
    if(itr == vector.end()) {
      std::cerr << "Couldn't find " << search << std::endl;
      return;
    }
    vector.erase(itr);
  }
}

#endif
