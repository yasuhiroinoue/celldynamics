/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "_class.h"

namespace output {
  void outputVTK(const Global*, unsigned int);
  void init_outputReconnection(const Global*);
  void outputReconnection(const Global*);
  void output_vertex(Global*);
}

#endif
