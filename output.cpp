/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "output.h"
#include "_parameters.h"
#include "_class.h"
#include <iostream>
#include <fstream>

namespace output {

void outputVTK(const Global *p_g, unsigned int step) {
  {
    char fname[100];
    sprintf(fname, "2dv_line%010u.vtk", step);
    std::ofstream fout(fname);

    //vtkファイルのヘッダ
    fout << "# vtk DataFile Version 2.0" << std::endl;
    fout << "2D-vertex" << std::endl;
    fout << "ASCII" << std::endl;
    fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

    //点の位置を書き込む
    fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
    for (int i = 0; i < (int)p_g->p_v.size(); i++) {
      Vertex *vp = p_g->p_v[i];
      fout << vp->loc[0].x << " ";
      fout << vp->loc[0].y << " ";
      fout << vp->loc[0].z << std::endl;
    }
    //線の要素を書き込む
    fout << "CELLS " << p_g->p_l.size() << " " << p_g->p_l.size() * 3 << std::endl;
    for (int i = 0; i < (int)p_g->p_l.size(); i++) {
      Line *lp = p_g->p_l[i];
      fout << "2 " << lp->vi[0] << " " << lp->vi[1] << std::endl;
    }
    //CELL_TYPES
    fout << "CELL_TYPES " << p_g->p_l.size() << std::endl;
    for (int i = 0; i < (int)p_g->p_l.size(); i++) {
      fout << "3" << std::endl;
    }

    fout << "CELL_DATA " <<  p_g->p_l.size() << std::endl;
    fout << "SCALARS line_tension float" << std::endl;
    fout << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < (int)p_g->p_l.size(); ++i) {
      auto *lp = p_g->p_l[i];
      fout << lp->lt << std::endl;
    }
  }
  {
    char fname[100];
    sprintf(fname, "2dv_face%010u.vtk", step);
    std::ofstream fout(fname);

    //vtkファイルのヘッダ
    fout << "# vtk DataFile Version 2.0" << std::endl;
    fout << "2D-vertex" << std::endl;
    fout << "ASCII" << std::endl;
    fout << "DATASET UNSTRUCTURED_GRID" << std::endl;

    //点の位置を書き込む
    fout << "POINTS " << p_g->p_v.size() << " float" << std::endl;
    for (int i = 0; i < (int)p_g->p_v.size(); i++) {
      Vertex *vp = p_g->p_v[i];
      fout << vp->loc[0].x << " ";
      fout << vp->loc[0].y << " ";
      fout << vp->loc[0].z << std::endl;
    }
    int cells_size = 0;

    for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
      Cellula *cp = p_g->p_c[i];
      cells_size += cp->vi.size() + 1;
    }

    fout << "CELLS " << p_g->p_c.size() << " " << cells_size << std::endl;
    for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
      Cellula *cp = p_g->p_c[i];
      fout << cp->vi.size();
      for (int j = 0; j < (int)cp->vi.size(); j++) {
        fout << " " << cp->vi[j];
      }
      fout << std::endl;
    }
    //CELL_TYPES
    fout << "CELL_TYPES " << p_g->p_c.size() << std::endl;
    for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
      fout << "7" << std::endl;
    }

    fout << "CELL_DATA " <<  p_g->p_c.size() << std::endl;
    fout << "SCALARS phase float" << std::endl;
    fout << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < (int)p_g->p_c.size(); i++) {
      auto *cp = p_g->p_c[i];
      double sin_t = sin(2.0 * M_PI * (cp->cell_time) / cp->cell_T - cp->cell_phase);
      fout << sin_t << std::endl;
    }
  }
}


void init_outputReconnection(const Global *p_g) {
  char fname[100];
  sprintf(fname, "2dv_reconnection.dat");
  std::ofstream fout(fname);
  fout << p_g->step << "	" << p_g->reconnection << std::endl;
}

void outputReconnection(const Global *p_g) {
  char fname[100];
  sprintf(fname, "2dv_reconnection.dat");
  std::ofstream fout(fname, std::ios::app);
  fout << p_g->step << "	" << p_g->reconnection << std::endl;
}

void output_vertex(Global *p_g) {
  std::ofstream ofs("vertex.dat");

  for(int i = 0; i < (int)p_g->p_v.size(); ++i) {
    ofs << p_g->p_v[i]->loc[0].x << " " << p_g->p_v[i]->loc[0].y << std::endl;
  }
}

}//namespace output
