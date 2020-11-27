/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <time.h>
#include <random>
#include <omp.h>

#include "vec.h"
#include "_class_and_variables.h"
#include "_class.h"
#include "_parameters.h"
#include "output.h"
#include "force.h"
#include "restructure.h"
#include "ODE_solver.h"
//#include "mymath.h"

void initPlain(Global*);
void calcCenter(Global*);
void yuragi(Global*);

//Global *p_g;

void readParameter(void) {
  // read "parameter.2dv" and read parameters for PCP
  char fna[100] = "parameter.2dv";

  std::ifstream fin(fna, std::ios::in);
  if (!fin.is_open()) {
    std::cout << "Error: missing parameter.2dv" << std::endl;
    exit(1);
  }

  fin >> K1_PCP_LENGTH;
  fin >> FLUCT;
  fin >> power_pcp;

  fin >> Pulse_T;
  fin >> phase_x;
  fin >> phase_y;
  unsigned long int _seed = 20110412;
  fin >> _seed;

  std::cout << "Line tension by PCP = " << K1_PCP_LENGTH << std::endl;
  std::cout << "Fluctuation = " << FLUCT << std::endl;
  std::cout << "Power of inner product (PCP) = " << power_pcp << std::endl;
  std::cout << "Pulse Period = " << Pulse_T << std::endl;
  std::cout << "Phase Shift along x-axis = " << phase_x << std::endl;
  std::cout << "Phase Shift along y-axis = " << phase_y << std::endl;

  Init_Random(_seed);
}


void initPlain(Global *p_g) {
  //頂点の座標を入れる。
  double y_tmp = 0.0;

  for (int j = 0; j < 2 * NUM_Y + 1; j++) {
    if (j % 2 == 0) {
      for (int i = 0; i < 2 * NUM_X; i++) {
        Vertex *tmp = new Vertex;
        if (i == 0) {
          tmp->loc[0].x = 0.5;
        }
        else if (i % 2 == 1) {
          tmp->loc[0].x = p_g->p_v.back()->loc[0].x + 1.0;
        }
        else if (i % 2 == 0) {
          tmp->loc[0].x = p_g->p_v.back()->loc[0].x + 2.0; //0.5 + (double)i*1.0;
        }
        else {
          std::cout << "Bug.Init" << std::endl;
          exit(0);
        }
        tmp->loc[0].y = y_tmp;
        tmp->loc[0].z = 0.0;

        p_g->p_v.push_back(tmp);
        //std::cout << i << " " << p_g->p_v[i]->loc.x << " " << p_g->p_v[i]->loc.y << std::endl;
      }
    }
    else if (j % 2 == 1) {
      for (int i = 0; i < 2 * NUM_X; i++) {

        Vertex *tmp = new Vertex;
        //tmp->loc.x = (double)i*1.0;
        if (i == 0) {
          tmp->loc[0].x = 0.0;
        }
        else if (i % 2 == 1) {
          tmp->loc[0].x = p_g->p_v.back()->loc[0].x + 2.0;
        }
        else if (i % 2 == 0) {
          tmp->loc[0].x = p_g->p_v.back()->loc[0].x + 1.0; //0.5 + (double)i*1.0;
        }
        else {
          std::cout << "Bug.Init" << std::endl;
          exit(0);
        }
        tmp->loc[0].y = y_tmp;
        tmp->loc[0].z = 0.0;

        p_g->p_v.push_back(tmp);
        //std::cout << i << " " << p_g->p_v[i]->loc.x << " " << p_g->p_v[i]->loc.y << std::endl;
      }
    }
    else {
      std::cout << "Bug.Init" << std::endl;
      exit(0);
    }
    y_tmp += sqrt(3.0) / 2.0;
  }

  //線に点のインデックスを入れる。
  for (int i = 0; i < NUM_X; i++) {
    for (int j = 0; j < 2 * NUM_Y; j++) {
      Line *tmp = new Line;
      tmp->vi[0] = 2 * i + 2 * NUM_X * j;
      tmp->vi[1] = 2 * i + 2 * NUM_X * (j + 1);
      tmp->lt = 0.0;

      p_g->p_l.push_back(tmp);
    }
    for (int j = 0; j < NUM_Y + 1; j++) {
      Line *tmp = new Line;
      tmp->vi[0] = 2 * i + 2 * 2 * NUM_X * j;
      tmp->vi[1] = 1 + 2 * i + 2 * 2 * NUM_X * j;
      tmp->lt = 0.0;

      p_g->p_l.push_back(tmp);
    }
    for (int j = 0; j < 2 * NUM_Y; j++) {
      Line *tmp = new Line;
      tmp->vi[0] = 1 + 2 * i + 2 * NUM_X * j;
      tmp->vi[1] = 1 + 2 * i + 2 * NUM_X * (j + 1);
      tmp->lt = 0.0;

      p_g->p_l.push_back(tmp);
    }
    if (i == NUM_X - 1) {
      break;
    }
    for (int j = 0; j < NUM_Y; j++) {
      Line *tmp = new Line;
      tmp->vi[0] = 2 * (NUM_X + 1) + 2 * i + 2 * 2 * NUM_X * j - 1;
      tmp->vi[1] = 2 * (NUM_X + 1) + 2 * i + 2 * 2 * NUM_X * j;
      tmp->lt = 0.0;

      p_g->p_l.push_back(tmp);
    }
  }

  //点に線のインデックスを入れる。
  for (int i = 0; i < (int)p_g->p_l.size(); i++) {
    Vertex *vp[2];
    vp[0] = p_g->p_v[p_g->p_l[i]->vi[0]];
    vp[1] = p_g->p_v[p_g->p_l[i]->vi[1]];
    vp[0]->li.push_back(i);
    vp[1]->li.push_back(i);
  }

  //細胞に線と点のインデックスを入れる。
  for (int i = 0; i < 2 * NUM_X - 1; i++) {
    if (i % 2 == 0) {
      for (int j = 0; j < NUM_Y; j++) {
        Cellula *tmp = new Cellula;

        //線のインデックスを入れる。
        tmp->li.push_back(2 * j + (6 * NUM_Y + 1)*i / 2);
        tmp->li.push_back(2 * j + 1 + (6 * NUM_Y + 1)*i / 2);
        tmp->li.push_back(2 * NUM_Y + j + (6 * NUM_Y + 1)*i / 2);
        tmp->li.push_back(2 * NUM_Y + (j + 1) + (6 * NUM_Y + 1)*i / 2);
        tmp->li.push_back((3 * NUM_Y + 1) + 2 * j + (6 * NUM_Y + 1)*i / 2);
        tmp->li.push_back((3 * NUM_Y + 1) + 2 * j + 1 + (6 * NUM_Y + 1)*i / 2);

        //点のインデックスを入れる。反時計回りに格納されるようにしている。
        tmp->vi.push_back(4 * NUM_X * j + i);
        tmp->vi.push_back(4 * NUM_X * j + i + 1);
        tmp->vi.push_back(2 * NUM_X + 4 * NUM_X * j + i + 1);
        tmp->vi.push_back(4 * NUM_X * (j + 1) + i + 1);
        tmp->vi.push_back(4 * NUM_X * (j + 1) + i);
        tmp->vi.push_back(2 * NUM_X + 4 * NUM_X * j + i);

        p_g->p_c.push_back(tmp);
      }
    }
    else if (i % 2 == 1) {
      for (int j = 0; j < NUM_Y - 1; j++) {
        Cellula *tmp = new Cellula;

        //線のインデックスを入れる。
        tmp->li.push_back(2 * j + (3 * NUM_Y + 1) + 1 + (6 * NUM_Y + 1) * (i - 1) / 2);
        tmp->li.push_back(2 * j + (3 * NUM_Y + 1) + 2 + (6 * NUM_Y + 1) * (i - 1) / 2);
        tmp->li.push_back(j + (5 * NUM_Y + 1) + (6 * NUM_Y + 1) * (i - 1) / 2);
        tmp->li.push_back((j + 1) + (5 * NUM_Y + 1) + (6 * NUM_Y + 1) * (i - 1) / 2);
        tmp->li.push_back(2 * j + 1 + (6 * NUM_Y + 1) * (i + 1) / 2);
        tmp->li.push_back(2 * j + 2 + (6 * NUM_Y + 1) * (i + 1) / 2);

        //点のインデックスを入れる。反時計回りに格納されるようにしている。
        tmp->vi.push_back(2 * NUM_X + 1 + 4 * NUM_X * j + (i - 1));
        tmp->vi.push_back(2 * NUM_X + 2 + 4 * NUM_X * j + (i - 1));
        tmp->vi.push_back(4 * NUM_X + 2 + 4 * NUM_X * j + (i - 1));
        tmp->vi.push_back(6 * NUM_X + 2 + 4 * NUM_X * j + (i - 1));
        tmp->vi.push_back(6 * NUM_X + 1 + 4 * NUM_X * j + (i - 1));
        tmp->vi.push_back(4 * NUM_X + 1 + 4 * NUM_X * j + (i - 1));

        p_g->p_c.push_back(tmp);
      }
    }
    else {
      std::cout << "Bug.Init" << std::endl;
      exit(0);
    }
  }

  //線に細胞のインデックスを入れる。
  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    Cellula *cp = p_g->p_c[i];
    for (int j = 0; j < (int)cp->li.size(); j++) {
      Line *lp = p_g->p_l[cp->li[j]];
      lp->ci.push_back(i);
    }
  }

  //点に細胞のインデックスを入れる。
  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    Cellula *cp = p_g->p_c[i];
    for (int j = 0; j < (int)cp->vi.size(); j++) {
      Vertex *vp = p_g->p_v[cp->vi[j]];
      vp->ci.push_back(i);
    }
  }

  for(int i = 0; i < (int)p_g->p_l.size(); ++i) {
    p_g->p_l[i]->K1_LENGTH = K1_LENGTH;
    p_g->p_l[i]->K1_PCP_LENGTH = K1_PCP_LENGTH;
    p_g->p_l[i]->K2_LENGTH = K2_LENGTH;
    p_g->p_l[i]->LENGTH_EQ = LENGTH_EQ;
  }

  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    p_g->p_c[i]->AREA_EQ = AREA_EQ;
    p_g->p_c[i]->K_AREA = K_AREA;
  }

  calcCenter(p_g);

  // ----- ad hoc -------------
  // 最下段の細胞と、最右端の細胞の、登録情報がおかしい。たぶん線に登録されている細胞数が実際より多い
  // とりあえず、動かない、リコネクションしないように固定するフラグをたてる

  //for( auto cp = begin( p_g->p_c); cp != end( p_g->p_c ); ++cp ){
  //	(*cp)->fix = 0;
  //if( (*cp)->center.y < 1.0 ) (*cp)->fix = 1;

  //if( (*cp)->center.x > 27.5 ) (*cp)->fix = 1;
  //if( (*cp)->center.y > 32.5 ) (*cp)->fix = 1;
  //}


  for ( auto &cp : p_g->p_c) {
    cp->fix = 0;
  }

  // -----------------------
  // ----- phase initialize ----
  int flag = 0;
  int flag_L = 0;
  int flag_S = 0;
  int interval = NUM_Y;

  double phase_tmp = 0.0;
  double phase_x_tmp = 0.0;


  for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {

    int count = (i + 1) - flag_L * NUM_Y - flag_S * (NUM_Y - 1);

    if ( count % interval == 0 ) {
      if ( flag == 0 ) {
        flag = 1;
        flag_L++;
        interval = NUM_Y - 1;
        phase_x_tmp += phase_x;
        phase_tmp = phase_x_tmp;
      }
      else if ( flag == 1 ) {
        flag = 0;
        flag_S++;
        interval = NUM_Y;
        phase_x_tmp += phase_x;
        phase_x_tmp = phase_x_tmp;
      }
    }

    Cellula *cp;
    cp = p_g->p_c[i];
    cp->cell_time = 0.0;
    cp->cell_phase = phase_tmp;
    cp->cell_T = Pulse_T;

    phase_tmp += phase_y;
  }


  //細胞に点のインデックスを入れる。
  /*
  for(int i=0; i<p_g->p_c.size(); i++){
  	Cellula *cp;
  	cp = p_g->p_c[i];
  	for(int j=0; j<cp->li.size(); j++){
  		Line *lp;
  		lp = p_g->p_l[cp->li[j]];
  		cp->vi.push_back(lp->vi[0]);
  		cp->vi.push_back(lp->vi[1]);
  	}
  	//重複の削除
  	sort(cp->vi.begin(), cp->vi.end());
  	cp->vi.erase(unique(cp->vi.begin(), cp->vi.end()), cp->vi.end());
  }*/
  //restruction::sortAntiClockwise();


  //for debug
	/*
  std::cout << "points" << std::endl;
	for(int i=0; i<(int)p_g->p_v.size(); i++){
		std::cout << i << " " << p_g->p_v[i]->loc[0].x << " " << p_g->p_v[i]->loc[0].y << std::endl;
		std::cout << "li ";
		for(int j=0; j<p_g->p_v[i]->li.size(); j++){
			std::cout << p_g->p_v[i]->li[j] << " ";
			if(j == p_g->p_v[i]->li.size()-1){
				std::cout << std::endl;
			}
		}
	}
	std::cout << "lines" << std::endl;
	for(int i=0; i<p_g->p_l.size(); i++){
		std::cout << i << " vi " << p_g->p_l[i]->vi[0] << " " << p_g->p_l[i]->vi[1] << std::endl;
		for(int j=0; j<p_g->p_l[i]->ci.size(); j++){
			std::cout << i << " ci " << p_g->p_l[i]->ci[j] << " ";
			if(j == p_g->p_l[i]->ci.size()-1){
				std::cout << std::endl;
			}
		}
	}
	std::cout << "cells" << std::endl;
	for(int i=0; i<p_g->p_c.size(); i++){
		std::cout << i << " li ";
		for(int j=0; j<p_g->p_c[i]->li.size(); j++){
			std::cout << p_g->p_c[i]->li[j] << " ";
			if(j == p_g->p_c[i]->li.size()-1){
				std::cout << std::endl;
			}
		}
		std::cout << i << " vi ";
		for(int j=0; j<p_g->p_c[i]->vi.size(); j++){
			std::cout << p_g->p_c[i]->vi[j] << " ";
			if(j == p_g->p_c[i]->vi.size()-1){
				std::cout << std::endl;
			}
		}
	}
  */
  /*
  for(int i = 0; i < p_g->p_c.size(); i++){
  	std::cout << i << " " << p_g->p_c[i]->center.x << " " << p_g->p_c[i]->center.y << std::endl;
  }*/
};

void calcCenter(Global *p_g) {
  //細胞に点が反時計回りで格納されていることを前提にしている。
  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    Cellula *cp = p_g->p_c[i];
    _vec<double> tmp = _vec<double>(0.0, 0.0, 0.0);
    double area_all = 0.0;
    for (int j = 0; j < (int)cp->vi.size(); j++) {
      double area_tmp;
      _vec<double> center_tmp;
      Vertex *vp[2];
      vp[0] = p_g->p_v[cp->vi[j]];
      vp[1] = p_g->p_v[cp->vi[(j + 1) % cp->vi.size()]];
      /*
      if (j != (int)cp->vi.size() - 1) {
        vp[1] = p_g->p_v[cp->vi[j + 1]];
      }
      else if (j == (int)cp->vi.size() - 1) { //bug fix j = cp --> j == cpに変更
        vp[1] = p_g->p_v[cp->vi[0]];
      }
      else {
        std::cout << "Bug.Cent" << std::endl;
        exit(0);
      }
      */

      area_tmp = 0.5 * (vp[0]->loc[0] % vp[1]->loc[0]).z;
      //area_tmp = 0.5 * (vp[0]->loc[0].x * vp[1]->loc[0].y - vp[1]->loc[0].x * vp[0]->loc[0].y);
      area_all += area_tmp;
      center_tmp = (vp[0]->loc[0] + vp[1]->loc[0]) / 3.0;
      tmp += area_tmp * center_tmp;
    }
    if ( area_all > 0.0 )
      cp->center = tmp / area_all;
    else {
      //cp->center /= tmp / (double)cp->vi.size();
      cp->center = tmp / area_all;
      std::cout << "cell area is less than 0: " << i << " " << cp->center.x << " " << cp->center.y << std::endl;
    }
  }
};
void updatePulse(Global *p_g) {
  //increment cell time by DELTA_TIME
  #pragma omp parallel for num_threads(THREAD_NUM)
  for ( int i = 0; i < (int)p_g->p_c.size(); i++ ) {
    Cellula *cp;
    cp = p_g->p_c[i];
    cp->cell_time += DELTA_TIME;
  }

};
void yuragi(Global *p_g) {
  /*srand((unsigned)time(NULL));
  Cellula *cp = p_g->p_c[rand()%(p_g->p_c.size())];
  for(int i = 0; i < cp->vi.size(); i++){
  	Vertex *vp = p_g->p_v[cp->vi[i]];
  	vp->loc[0].x -= 0.2;
  }*/

  //inoue added 2016.10.25
  for (int i = 0; i < (int)p_g->p_v.size(); i++) {
    Vertex *vp = p_g->p_v[i];
    _vec<double> xi = _vec<double>(RAND() - 0.50, RAND() - 0.50, RAND() - 0.50);
    xi *= 2.0 * sqrt(3.0) * FLUCT * sqrt(DELTA_TIME);
    vp->loc[0] += xi;
  }


  //
};

void lt_reset_for_viz(Global *p_g) {
  //reset line tension into zero
  for (int i = 0; i < (int)p_g->p_l.size(); i++) {
    Line *lp = p_g->p_l[i];
    lp->lt = 0.0;
  }
}

void topologyCheck(const Global *p_g) {
  // if a vertex has more than 3 edges, the program stops
  int flag = 0;
  for ( int i = 0; i < (int)p_g->p_v.size(); i++ ) {
    Vertex *vp = p_g->p_v[i];
    if ( vp->ci.size() > 3 ) flag++;
  }
  if ( flag > 0 ) {
    std::cout << "more than three edges" << std::endl;
    exit(0);
  }


  for ( auto cp = begin(p_g->p_c); cp != end(p_g->p_c); ++cp ) {
    (*cp)->flag_rec = 0;
  }
}

/*
double totalEnergy(const Global *p_g) {
  double total_energy = 0;
  for(int cidx = 0; cidx < (int)(p_g->p_c.size()); ++cidx) {
    double circumference = 0;
    for(int lidx : p_g->p_c[cidx]->li) {
      circumference += (p_g->p_l[lidx]->vi[0] - p_g->p_l[lidx]^>vi[1]).norm();
    }
  }
}
*/

int main(void) {
  Global *p_g = new Global;
  readParameter();
  initPlain(p_g);
  //yuragi(p_g);
  output::outputVTK(p_g, 0);

  output::init_outputReconnection(p_g);

  int STEP_RECONNECT = (int)(TIME_RECONNECT / DELTA_TIME);

  int division_step = TIME_CELL_DIVISION / DELTA_TIME;

  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_int_distribution<> rand_cidx(0, p_g->p_c.size() - 1);
  std::uniform_real_distribution<> rand_axis(0, 2.0 * M_PI);
  std::uniform_real_distribution<> cell_time_init(0, division_step * DELTA_TIME);

  for(int cidx = 0; cidx < (int)p_g->p_c.size(); ++cidx) {
    p_g->p_c[cidx]->cell_time = cell_time_init(mt);
  }

  for (unsigned int num = 1; num <= STEP_END; num++) {
    p_g->step++;

    //	std::cout << "center" << std::endl;
    calcCenter(p_g);

    //
    if(num % STEP_RECONNECT == (unsigned int)STEP_RECONNECT / 2) {
      restructure::cellIntersection(p_g);
      calcCenter(p_g);
      restructure::removeTriangleVoid(p_g);
    }
    //

    //	std::cout << "rearrange" << std::endl;
    if ( (p_g->step % STEP_RECONNECT) == 0 ) {
      restructure::cellRearrange2(p_g);
      calcCenter(p_g);
    }

	/*
    if(num == 1000 ){//% 1000001 == 0) {
      int apoptosis_cidx = rand_cidx(mt);
      std::cout << "apoptosis cidx-" << apoptosis_cidx << std::endl;
      for(int lidx : p_g->p_c[apoptosis_cidx]->li) {
        p_g->p_l[lidx]->K2_LENGTH *= 1.5;
        //p_g->p_l[lidx]->K1_LENGTH *= 0.90;
      }
    }
	*/
    /*
    if(num % 1000001 == 0) {
      int apoptosis_cidx = rand_cidx(mt);
      std::cout << "apoptosis cidx-" << apoptosis_cidx << std::endl;
      p_g->shrinking_ci.push_back(apoptosis_cidx);
    }
    if(num % 103 == 0) {
      restructure::apoptosis(p_g, 0.99);
    }
    */


    int cidx = rand_cidx(mt);
    _vec<double> axis;
    double theta = rand_axis(mt);
    axis.x = cos(theta);
    axis.y = sin(theta);
    axis.z = 0.0;
    if(p_g->p_c[cidx]->cell_time / DELTA_TIME >= division_step && restructure::isConvex(p_g, cidx)) {
      try {
        restructure::cellDivision(p_g, cidx, axis);
        p_g->p_c[cidx]->cell_time = 0.0;
        calcCenter(p_g);
        std::cout << "  divide cellula " << cidx << std::endl;
      } catch (char const *errmsg){
        //std::cout << errmsg << std::endl;
      }
    }
    
    topologyCheck(p_g);

    //	std::cout << "viz reset" << std::endl;
    lt_reset_for_viz(p_g);

    //	std::cout << "update pulse" << std::endl;
    updatePulse(p_g);

    //	std::cout << "motion" << std::endl;
    if (DEGREE_ACCURACY == 1) {
      ODE_solver::motionVertexFirst(p_g);
      //yuragi(& p_g);
    }
    else if (DEGREE_ACCURACY == 2) {
      ODE_solver::motionVertexSecond(p_g);
    }
    else {
      std::cout << "Error: DEGREE_ACCURACY must be 1 or 2." << std::endl;
      exit(0);
    }

    if (p_g->step % PERIOD_PARAVIEW == 0) {
      //		std::cout << "vtk" << std::endl;
      (void)output::outputVTK(p_g, p_g->step);
      output::outputReconnection(p_g);
      //		std::cout << "vtk end" << std::endl;
    }
  }

  std::cout << "finished!" << std::endl;
  return 0;
}
