/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include <cmath>
#include "_parameters.h"

/*****  parameters  *****************************************************/
//物理パラメータ ※粘性係数ηは１とする。
// ---- 固定したパラメータ ---------- 細胞六角形　一辺の長さ1 になるように決めた
double K_AREA = 1e+0;			//面積に対する弾性係数
double K2_LENGTH = 1.0;			//周長に対する弾性係数
double K1_LENGTH = -6.0;		//周長に対する線張力 //     -6.0で辺の長さ1で釣り合い   //-4.0(smaller hexagon)//-8.0(soft)
double LENGTH_EQ = 0.0;//6.0		//今回は平衡周長が定数のためマクロにしておく
double AREA_EQ = 2.60;			//今回は平衡面積が定数のためマクロにしておく
// --------------------- PCPの部分 ----------------------------
double	K1_PCP_LENGTH =	0.2;		//PCP軸に沿った線張力 0 < K1_PCP_LENGTH < |K1_LENGTH| 0以上でK1_LENGTHの絶対値より小さい
double	FLUCT  = 0e+0;			//ゆらぎの大きさ
double	power_pcp = 2.0;		//PCP軸との内積(cos(theta))をpower_pcp乗する。

double	Pulse_T = 1.0;			//細胞のパルス
double	phase_y = (0.25 * M_PI);		//初期位相 y方向のズレ
double	phase_x = (0.50 * M_PI);		//初期位相 x方向のずれ
//------------------------------------------------------------
/************************************************************************/
