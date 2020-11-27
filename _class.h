/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _CLASS_H
#define _CLASS_H

#include <iostream>
#include <vector>
#include "_parameters.h"
#include "vec.h"

class Vertex {
 public:
  //int vid;
  std::vector<int> li;		//辺のインデックス
  std::vector<int> ci;		//細胞のインデックス
  _vec<double> loc[DEGREE_ACCURACY];	//位置ベクトル
  _vec<double> frc[DEGREE_ACCURACY];	//力のベクトル
  _vec<double> frc_thread[THREAD_NUM];	//並列計算用に、各コアで力を代入
  void debug_print();
};

class Line {
 public:
  int vi[2];		//頂点のインデックス
  std::vector<int> ci;		//細胞のインデックス
  double lt {};		//line tension
  double lt_thread[THREAD_NUM];
  double K1_LENGTH;
  double K1_PCP_LENGTH;
  double K2_LENGTH;
  double LENGTH_EQ;
  void debug_print();
};

class Cellula {
 public:
  std::vector<int> vi;	//頂点のインデックス
  std::vector<int> li;	//線のインデックス
  _vec<double> center;	//重心座標

  double cell_time {};	//細胞時間
  double cell_phase {};	//細胞時間の初期位相
  double cell_T {};	//細胞のパルス周期

  int flag_rec {};	//そのステップでリコネクションしたか//同一ステップで2回以上させない
  int fix {};		//細胞の固定
  double AREA_EQ; //面積
  double K_AREA;

  void debug_print();
};

class Global {
  // システム全体を表すクラス
  public:
  unsigned int step;
  std::vector<Cellula *> p_c;
  std::vector<Line *> p_l;
  std::vector<Vertex *> p_v;

  std::vector<int> shrinking_ci;

  unsigned int reconnection; //the number of times to reconnect

  Global();
};

#endif
