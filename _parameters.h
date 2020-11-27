/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _PRM_H
#define _PRM_H
/*****  parameters  *****************************************************/
//並列計算パラメータ
constexpr int THREAD_NUM = 4;		//使用するコア数
//初期形状パラメータ
constexpr int NUM_X = 5;			//初期形状を決定。x方向
constexpr int NUM_Y = 5;  			//初期形状を決定。y方向

//システムパラメータ
constexpr double DELTA_TIME = 1e-04;	//1ステップに相当する時間間隔
constexpr int PERIOD_PARAVIEW = 10000;	//paraview出力間隔
constexpr int STEP_END = 10000000;		//計算上限
constexpr int DEGREE_ACCURACY = 2;		//計算精度。2次精度まで。2次にすると計算時間が長くなる。

//物理パラメータ ※粘性係数ηは１とする。
// ---- 固定したパラメータ ---------- 細胞六角形　一辺の長さ1 になるように決めた
extern double K_AREA;			//面積に対する弾性係数
extern double K2_LENGTH;			//周長に対する弾性係数
extern double K1_LENGTH;		//周長に対する線張力 //     -6.0で辺の長さ1で釣り合い   //-4.0(smaller hexagon)//-8.0(soft)
extern double LENGTH_EQ;//6.0		//今回は平衡周長が定数のためマクロにしておく
extern double AREA_EQ;			//今回は平衡面積が定数のためマクロにしておく

//reconnection用パラメータ
constexpr double L_THRESHOLD = 1e-01;	//reconnectionの閾値
constexpr double L_RECONNECTED = 5e-02;	//reconnection後、長さは閾値＋この値となる。
constexpr double TIME_RECONNECT = 10.0;	//reconnectionのタイミング
// -----------------------------------

//Cell division
constexpr double TIME_CELL_DIVISION = 1000.0;//cell cycle duration
//------------------------------------

// --------------------- PCPの部分 ----------------------------
extern double	K1_PCP_LENGTH;		//PCP軸に沿った線張力 0 < K1_PCP_LENGTH < |K1_LENGTH| 0以上でK1_LENGTHの絶対値より小さい
extern double	FLUCT;			//ゆらぎの大きさ
extern double	power_pcp;		//PCP軸との内積(cos(theta))をpower_pcp乗する。

extern double	Pulse_T;			//細胞のパルス
extern double	phase_y;		//初期位相 y方向のズレ
extern double	phase_x;		//初期位相 x方向のずれ
//------------------------------------------------------------
/************************************************************************/

#endif
