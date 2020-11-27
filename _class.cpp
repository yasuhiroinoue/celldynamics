/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "_class.h"
#include <iostream>
#include <vector>

Global::Global() {
  Global::step = 0;
  Global::reconnection = 0;
}

void Vertex::debug_print() {
  std::cerr << "lidx: ";
  for(int lidx: Vertex::li) {
    std::cerr << " " << lidx;
  }
  std::cerr << std::endl;

  std::cerr << "cidx: ";
  for(int cidx: Vertex::ci) {
    std::cerr << " " << cidx;
  }
  std::cerr << std::endl;
}

void Line::debug_print() {
  std::cerr << "vidx: " << Line::vi[0] << " " << Line::vi[1] << std::endl;

  std::cerr << "cidx: ";
  for(int cidx: Line::ci) {
    std::cerr << " " << cidx;
  }
  std::cerr << std::endl;
}

void Cellula::debug_print() {
  std::cerr << "vidx:";
  for(int vidx: Cellula::vi) {
    std::cerr << " " << vidx;
  }
  std::cerr << std::endl;

  std::cerr << "lidx: ";
  for(int lidx: Cellula::li) {
    std::cerr << " " << lidx;
  }
  std::cerr << std::endl;

}
