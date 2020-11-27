/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "ODE_solver.h"
#include "_parameters.h"
#include "force.h"

namespace ODE_solver {
//1次精度で微分方程式を解く（オイラー法）。
void motionVertexFirst(Global *p_g) {
//
  force::calcLineForce(p_g, 0);
  force::calcAreaForce(p_g, 0);


  force::OMP_Reduction_Frc(p_g, 0, 1);
  force::OMP_Reduction_Lt(p_g);


  #pragma omp parallel for num_threads(THREAD_NUM)
  for (int i = 0; i < (int)p_g->p_v.size(); i++) {
    Vertex *vp = p_g->p_v[i];

    int fix_flag = 0;
    for ( int j = 0; j < (int)vp->ci.size(); j++ ) {
      fix_flag += p_g->p_c[vp->ci[j]]->fix;
    }
    if ( fix_flag == 0 )
      vp->loc[0] += vp->frc[0] * DELTA_TIME;

    vp->frc[0] = _vec<double>(0.0, 0.0, 0.0);
  }
};

//2次精度で微分方程式を解く(中点法)。ただし、中点の座標と力は1次精度で求める。
void motionVertexSecond(Global *p_g) {
  //Δt/2後の位置を求める。
  force::calcLineForce(p_g, 0);
  force::calcAreaForce(p_g, 0);

  force::OMP_Reduction_Frc(p_g, 0, 2);

  #pragma omp parallel for num_threads(THREAD_NUM)
  for (int i = 0; i < (int)p_g->p_v.size(); i++) {
    Vertex *vp = p_g->p_v[i];
    vp->loc[1] = vp->loc[0] + vp->frc[0] * DELTA_TIME / 2.0;
    vp->frc[0] = _vec<double>(0.0, 0.0, 0.0);
  }
  //Δt/2後の力を求め、修正したΔt後の位置を出す。
  force::calcLineForce(p_g, 1);
  force::calcAreaForce(p_g, 1);


  force::OMP_Reduction_Frc(p_g, 1, 2);
  force::OMP_Reduction_Lt(p_g);

  #pragma omp parallel for num_threads(THREAD_NUM)
  for (int i = 0; i < (int)p_g->p_v.size(); i++) {
    Vertex *vp = p_g->p_v[i];
    vp->loc[0] += vp->frc[1] * DELTA_TIME;
    vp->frc[1] = _vec<double>(0.0, 0.0, 0.0);
  }
};

} // namespce ODE_solver
