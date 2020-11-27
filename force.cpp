/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "force.h"
#include "_class.h"
#include "_parameters.h"
#include "_class_and_variables.h"
#include "vec.h"
#include <omp.h>

namespace force {
void calcLineForce(Global *p_g, int deg) {
  //deg: degree accuracy
  #pragma omp parallel for num_threads(THREAD_NUM)
  for (int i = 0; i < (int)p_g->p_c.size(); i++) {

    int id_tnum = (int)omp_get_thread_num();

    double i_length = 0.0;
    Cellula *cp = p_g->p_c[i];
    //i番目の細胞の周長を計算する。
    for (int j = 0; j < (int)cp->li.size(); j++) {
      Line *lp = p_g->p_l[cp->li[j]];
      Vertex *vp[2];
      vp[0] = p_g->p_v[lp->vi[0]];
      vp[1] = p_g->p_v[lp->vi[1]];

      _vec<double> ij_line = vp[1]->loc[deg] - vp[0]->loc[deg];
      double ij_length = ij_line.norm();

      i_length += ij_length;
    }

    //線ごとに頂点を引っ張ってきて勾配と力を計算する。弾性力
    for (int j = 0; j < (int)cp->li.size(); j++) {
      Line *lp = p_g->p_l[cp->li[j]];
      Vertex *vp[2];
      vp[0] = p_g->p_v[lp->vi[0]];
      vp[1] = p_g->p_v[lp->vi[1]];

      _vec<double> edge = vp[0]->loc[deg] - vp[1]->loc[deg];
      double edge_length = edge.norm();
      _vec<double> l_grad;
      if ( edge_length > 1e-5 ) l_grad = edge / edge_length;
      _vec<double> frc_tmp = l_grad * lp->K2_LENGTH * (i_length - lp->LENGTH_EQ); //2次ポテンシャル(弾性エネルギー)

      vp[0]->frc_thread[id_tnum] -= frc_tmp;
      vp[1]->frc_thread[id_tnum] += frc_tmp;
    }

    //線張力

    for (int j = 0; j < (int)cp->li.size(); j++) {
      Line *lp = p_g->p_l[cp->li[j]];
      Vertex *vp[2];
      vp[0] = p_g->p_v[lp->vi[0]];
      vp[1] = p_g->p_v[lp->vi[1]];

      _vec<double> edge = vp[0]->loc[deg] - vp[1]->loc[deg];
      double edge_length = edge.norm();

      _vec<double> l_grad;
      if ( edge_length > 1e-5 ) l_grad = edge / edge_length;


      //1次ポテンシャル(線張力エネルギー) addedd by inoue 2016.10.25
      // Circumference line tension ------------------------
      _vec<double> frc_tmp = l_grad * lp->K1_LENGTH;

      // PCP regulation ------------------------------------
      _vec<double> pcp_vec = _vec<double>(1.0, 0.0, 0.0); //PCPによるx軸に沿った収縮
      pcp_vec /= pcp_vec.norm();//単位ベクトルにする

      double div = l_grad * pcp_vec; //div = cos(theta)のことです
      //frc_tmp += (div*div*div*div)*l_grad*K1_PCP_LENGTH;
      double sint_t = sin(2.0 * M_PI * (cp->cell_time) / cp->cell_T - cp->cell_phase);
      double line_tension_pcp = 0.0;

      //if( lp->ci.size() != 1 ){//境界型エッジには作用させない
      line_tension_pcp = pow(div, power_pcp) * lp->K1_PCP_LENGTH * (sint_t + 1.0);
      //}

      frc_tmp += line_tension_pcp * l_grad;

      // Fluctuation ---------------------------------------
      //_vec<double> xi = _vec<double>(RAND-0.50,RAND-0.50,0.0); xi *= 2.0*sqrt(3.0)*FLUCT;
      //frc_tmp += xi;


      // Paraview VTK 表示用 line tensionの大きさ ----------
      lp->lt_thread[id_tnum] += line_tension_pcp;

      //----------------- ここまで -------------------------


      vp[0]->frc_thread[id_tnum] -= frc_tmp;
      vp[1]->frc_thread[id_tnum] += frc_tmp;
    }

    /*ややこしいコードになった上に挙動がおかしいのでコメントアウト
    //i番目の細胞のうちj番目の頂点を引っ張ってくる。
    for(int j=0; j<cp->vi.size(); j++){
    	Vertex *vp[2];
    	vp[0] = p_g->p_v[cp->vi[j]];	//vp[0]がその点を表す。vp[1]は勾配の計算に使う。
    	_vec<double> l_grad = _vec<double>(0.0, 0.0, 0.0);
    	//j番目の点が構成する線についてループさせることにより、各点における勾配を計算
    	for(int k=0; k<vp[0]->li.size(); k++){
    		Line *lp = p_g->p_l[vp[0]->li[k]];

    		//i番目の細胞にない線を参照している場合は計算しない。
    		bool flag;
    		for(int l=0; l<lp->ci.size(); l++){
    			if(lp->ci[l] == i){
    				flag = true;
    				break;
    			}
    			if(l == lp->ci.size()-1){
    				flag = false;
    			}
    		}
    		if(flag == false){
    			break;
    		}

    		//基準点と逆側の点を引っ張ってくる。
    		if(p_g->p_v[lp->vi[0]] != vp[0]){
    			vp[1] = p_g->p_v[lp->vi[0]];
    		} else if(p_g->p_v[lp->vi[1]] != vp[0]){
    			vp[1] = p_g->p_v[lp->vi[1]];
    		} else{
    			std::cout << "Bug." << std::endl;
    			exit(0);
    		}

    		_vec<double> dist_vec = vp[0]->loc - vp[1]->loc;	//距離ベクトル
    		double dist_length = dist_vec.norm();

    		l_grad = dist_vec/dist_length;	//+=
    		vp[0]->frc += (-1.0)*l_grad*K2_LENGTH*(i_length-LENGTH_EQ);
    	}
    	//i番目の細胞内のポテンシャルによって点に加わる力を計算する。
    	//vp[0]->frc += (-1.0)*l_grad*K2_LENGTH*(i_length-LENGTH_EQ);
    	//std::cout << vp[0]->frc.x << std::endl;
    }
    */
  }
};

void calcAreaForce(Global *p_g, int deg) {
  //deg: degree accuracy
  #pragma omp parallel for num_threads(THREAD_NUM)
  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    int id_tnum = (int)omp_get_thread_num();

    Cellula *cp = p_g->p_c[i];
    double i_area = 0.0;
    //i番目の細胞の面積を計算。点が反時計回りに格納されていることを前提にしている。
    for (int j = 0; j < (int)cp->vi.size(); j++) {
      Vertex *vp[2];
      vp[0] = p_g->p_v[cp->vi[j]];
      if (j != (int)cp->vi.size() - 1) {
        vp[1] = p_g->p_v[cp->vi[j + 1]];
      }
      else if (j == (int)cp->vi.size() - 1) {
        vp[1] = p_g->p_v[cp->vi[0]];
      }
      else {
        std::cout << "Bug.Area" << std::endl;
        exit(0);
      }
      i_area += 0.5 * (vp[0]->loc[deg].x * vp[1]->loc[deg].y - vp[1]->loc[deg].x * vp[0]->loc[deg].y);
    }

	

    //面積勾配を計算。
    for (int j = 0; j < (int)cp->vi.size(); j++) {
      Vertex *vp[3];
      vp[1] = p_g->p_v[cp->vi[j]];

      if (j == 0) {
        vp[0] = p_g->p_v[cp->vi[cp->vi.size() - 1]];
      }
      else {
        vp[0] = p_g->p_v[cp->vi[j - 1]];
      }

      if (j == (int)cp->vi.size() - 1) {
        vp[2] = p_g->p_v[cp->vi[0]];
      }
      else {
        vp[2] = p_g->p_v[cp->vi[j + 1]];
      }

      _vec<double> s_grad = _vec<double>(0.0, 0.0, 0.0);
      s_grad.x = 0.5 * (vp[2]->loc[deg].y - vp[0]->loc[deg].y);
      s_grad.y = 0.5 * (vp[0]->loc[deg].x - vp[2]->loc[deg].x);

      _vec<double> frc_tmp = (-1.0) * cp->K_AREA * (i_area - cp->AREA_EQ) * s_grad;
      vp[1]->frc_thread[id_tnum] += frc_tmp;
    }
  }
};

void OMP_Reduction_Frc(Global *p_g, int deg, int acc) { //accは時間精度
  // deg: degree accuracy
  #pragma omp parallel for num_threads(THREAD_NUM)
  for ( int i = 0; i < (int)p_g->p_v.size(); i++ ) {
    Vertex *vp = p_g->p_v[i];
    vp->frc[deg] = _vec<double>(0., 0., 0.);
    for ( int tnum = 0; tnum < THREAD_NUM; tnum++ ) {
      vp->frc[deg] += vp->frc_thread[tnum];
      vp->frc_thread[tnum] = _vec<double>(0.0, 0.0, 0.0);
    }
    //Fluctuation --------------------------------------- １次精度、２次精度でちゃんとわける.２次で両方に入れると理論上おかしい
    if ( acc == 1 && deg == 0) {
      _vec<double> xi = _vec<double>(RAND() - 0.50, RAND() - 0.50, 0.0);
      xi *= 2.0 * sqrt(3.0) * FLUCT;
      vp->frc[deg] += xi;
    }
    else if ( acc == 2 && deg == 1 ) {
      _vec<double> xi = _vec<double>(RAND() - 0.50, RAND() - 0.50, 0.0);
      xi *= 2.0 * sqrt(3.0) * FLUCT;
      vp->frc[deg] += xi;
    }
  }
}

void OMP_Reduction_Lt(Global *p_g) {
  #pragma omp parallel for num_threads(THREAD_NUM)
  for ( int i = 0; i < (int)p_g->p_l.size(); i++ ) {
    Line *lp = p_g->p_l[i];
    lp->lt = 0.0;
    for ( int tnum = 0; tnum < THREAD_NUM; tnum++ ) {
      lp->lt += lp->lt_thread[tnum];
      lp->lt_thread[tnum] = 0.0;
    }
  }
}

} // namespace force
