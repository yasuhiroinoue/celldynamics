/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _FRC_H
#define _FRC_H

#include "_class.h"

namespace force {
  void calcLineForce(Global*, int);
  void calcAreaForce(Global*, int);
  void OMP_Reduction_Frc(Global*, int, int);
  void OMP_Reduction_Lt(Global*);
}

#endif
