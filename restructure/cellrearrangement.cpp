/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include "../restructure.h"
#include "../_class.h"
#include "../_parameters.h"
#include "../vec.h"
#include <iostream>
#include <vector>
#include <algorithm>

namespace restructure {

void sortAntiClockwise(Global *p_g) {
//	std::cout << "Anticlock start" << std::endl;

  for (int i = 0; i < (int)p_g->p_c.size(); i++) {
    Cellula *cp = p_g->p_c[i];
    std::vector<int> vi_tmp;

    //std::cout << "First section" << std::endl;
    //時計回りまたは反時計回りになるようにソートする。
    vi_tmp.push_back(p_g->p_l[cp->li[0]]->vi[0]);

    int li0_0 = p_g->p_l[cp->li[0]]->vi[0];
    int li0_1 = p_g->p_l[cp->li[0]]->vi[1];

    for (int j = 1; j < (int)cp->li.size();  j++) {
      for (int k = j; k < (int)cp->li.size(); k++) {
        Line *lp = p_g->p_l[cp->li[k]];

        if (lp->vi[0] == vi_tmp.back()) {

          vi_tmp.push_back(lp->vi[1]);
          int tmp = cp->li[j];
          cp->li[j] = cp->li[k];
          cp->li[k] = tmp;

          break;

        }
        else if (lp->vi[1] == vi_tmp.back()) {


          vi_tmp.push_back(lp->vi[0]);
          int tmp = cp->li[j];
          cp->li[j] = cp->li[k];
          cp->li[k] = tmp;
          break;

        }
        else {
          continue;
        }


      }
    }

    //重複判定
    /*
    int flag = 0;
    for( int i = 0; i < vi_tmp.size()-1; i++ ){
    	for( int j = (i+1); j < vi_tmp.size(); j++ ){
    		if( vi_tmp[i] == vi_tmp[j] ){
    			std::cout << vi_tmp[i] << " == " << vi_tmp[j] << std::endl;
    			flag++;
    		}
    	}
    }
    if( flag > 0 ){ std::cout << flag << std::endl; exit(0); }
    */


    //std::cout << "Second section " << cp->vi.size() << " " << vi_tmp.size() << std::endl;
    //符号付面積の符号で時計回りか反時計回りか判断できる。
    double area_tmp = 0.0;
    for (int j = 0; j < (int)vi_tmp.size(); j++) {
      Vertex *vp[2];
      vp[0] = p_g->p_v[vi_tmp[j]];
      if (j != (int)cp->vi.size() - 1) {
        vp[1] = p_g->p_v[vi_tmp[j + 1]];
      }
      else if (j == (int)cp->vi.size() - 1) {
        vp[1] = p_g->p_v[vi_tmp[0]];
      }
      else {
        std::cout << "Bug.Anticlock" << std::endl;
        exit(0);
      }
      area_tmp += 0.5 * (vp[0]->loc[0].x * vp[1]->loc[0].y - vp[1]->loc[0].x * vp[0]->loc[0].y);
    }


    //std::cout << "Third section" << std::endl;
    //正なら反時計回りなのでそのまま。負なら時計回りなので反転する。
    if (area_tmp > 0) {
      cp->vi.clear();
      for (int j = 0; j < (int)vi_tmp.size(); j++) {
        cp->vi.push_back(vi_tmp[j]);
      }
    }
    else if (area_tmp < 0) {
      reverse(begin(vi_tmp), end(vi_tmp));
      cp->vi.clear();
      for (int j = 0; j < (int)vi_tmp.size(); j++) {
        cp->vi.push_back(vi_tmp[j]);
      }
    }
    else {
      std::cout << "Error: area must not be 0." << std::endl;
      exit(0);
    }
  }

//	std::cout << "Anticlock end" << std::endl;
};

void cellRearrange(Global *p_g) {

  bool rc_flag = false;

  for (int i = 0; i < (int)p_g->p_l.size(); i++) {
    Line *lp = p_g->p_l[i];
    Vertex *vp[2];
    vp[0] = p_g->p_v[lp->vi[0]];
    vp[1] = p_g->p_v[lp->vi[1]];
    _vec<double> norm_vec = vp[0]->loc[0] - vp[1]->loc[0];
    double norm = norm_vec.norm();

    //線の長さが閾値より少ないとき、繋ぎ変え判定を行う
    if (norm < L_THRESHOLD ) {
      std::vector<int> cell_rec;

      //繋ぎ変え許可のトポロジー判定
      int pre_flag = 0;
      if ( lp->ci.size() == 1 && p_g->p_c[lp->ci[0]]->flag_rec == 0 &&  p_g->p_c[lp->ci[0]]->fix == 0) { //境界型のエッジ
        if ( vp[0]->ci.size() == 2 && vp[1]->ci.size() == 2 ) { //端点に、細胞が２個と２個
          pre_flag = 3;
        }
      }
      else if ( lp->ci.size() == 2 ) {  //共有型のエッジ
        if ( p_g->p_c[lp->ci[0]]->flag_rec == 0 && p_g->p_c[lp->ci[1]]->flag_rec == 0 &&  p_g->p_c[lp->ci[0]]->fix == 0 &&  p_g->p_c[lp->ci[1]]->fix == 0) {
          int prod = vp[0]->ci.size() * vp[1]->ci.size();
          int sum =  vp[0]->ci.size() + vp[1]->ci.size();
          if ( prod == 6 && sum == 5) pre_flag = 2; //端点に、細胞が２個と３個
          else if ( prod == 9 && sum == 6 ) { //端点に、細胞が３個と３個
            int post_c0 = 0, post_c1 = 0;

            for ( auto cc : vp[0]->ci ) {
              if ( cc != lp->ci[0] && cc != lp->ci[1] ) {
                post_c0 = cc;
                //std::cout << cc << " " << lp->ci[0] << " " << lp->ci[1] << std::endl;
              }
            }

            for ( auto cc : vp[1]->ci ) {
              if ( cc != lp->ci[0] && cc != lp->ci[1] ) {
                post_c1 = cc;
              }
            }

            _vec<double> r_c0c1 = p_g->p_c[post_c0]->center - p_g->p_c[post_c1]->center;
            _vec<double> r_v0v1 = vp[0]->loc[0] - vp[1]->loc[0];

            if ( r_c0c1 * r_v0v1 > 0 ) { //交差していない
              pre_flag = 1;
            }
          }
        }
      }

      if ( pre_flag == 1 ) {
        p_g->reconnection++;
        cell_rec.push_back(lp->ci[0]);
        cell_rec.push_back(lp->ci[1]);
        //std::cout << "pre_flag = 1 start" << std::endl;
        rc_flag = true;
        //座標の移動
        _vec<double> mid_vec = 0.5 * (vp[0]->loc[0] + vp[1]->loc[0]);
        _vec<double> tmp_v0 = _vec<double>(0.0, 0.0, 0.0);
        _vec<double> tmp_v1 = _vec<double>(0.0, 0.0, 0.0);
        tmp_v0.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[0]->loc[0].y);
        tmp_v0.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[0]->loc[0].x - mid_vec.x);
        tmp_v1.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[1]->loc[0].y);
        tmp_v1.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[1]->loc[0].x - mid_vec.x);

        vp[0]->loc[0] = mid_vec + tmp_v0;
        vp[1]->loc[0] = mid_vec + tmp_v1;

        //線が新しく持つ細胞のインデックス
        std::vector<int> ci_tmp;
        for (int j = 0; j < 2; j++) {
          for (auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++) {
            bool flag = false;
            for (auto jt = begin(lp->ci); jt != end(lp->ci); jt++) {
              if ((*it) == (*jt)) {
                break;
              }
              if (jt == end(lp->ci) - 1) {
                flag = true;
              }
            }
            if (flag == true) {
              ci_tmp.push_back(*it);
            }
          }
        }


        std::cout << "pre_flag = " << pre_flag << " cell size = " << ci_tmp.size() << std::endl;

        //細胞の重心との位置関係で移動後の点がどの位置にあるのか判断する。細胞の変形があまりにも激しいとバグるかもしれない。
        Cellula *cp0 = p_g->p_c[lp->ci[0]];
        Cellula *cp1 = p_g->p_c[lp->ci[1]];
        _vec<double> c_vec0 = cp0->center - mid_vec;
        //もし内積が正（重心と同じ側にvp[0]がある）ならばcp0がvp0の移動先
        if (c_vec0 * tmp_v0 > 0) {
          //点0側の細胞から点1の情報を消す。両方の点についてやる。
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[1]) {
              cp0->vi.erase(it);
              break;
            }
          }
          for (auto it = begin(cp1->vi); it != end(cp1->vi); it++) {
            if (*it == lp->vi[0]) {
              cp1->vi.erase(it);
              break;
            }
          }

        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[0]) {
              cp0->vi.erase(it);
              break;
            }
          }
          for (auto it = begin(cp1->vi); it != end(cp1->vi); it++) {
            if (*it == lp->vi[1]) {
              cp1->vi.erase(it);
              break;
            }
          }

        }
        else {
          std::cout << "Error: recconection was failed. c_vec0*tmp_v0 = " << c_vec0 *tmp_v0 << " " << cp0->center.x << std::endl;
          exit(0);
        }

        //細胞から線のインデックスを消去する。
        for (int j = 0; j < (int)lp->ci.size(); j++) {
          Cellula *cp = p_g->p_c[lp->ci[j]];
          for (auto it = begin(cp->li); it != end(cp->li); it++) {
            if (*it == i) {
              cp->li.erase(it);
              break;
            }
          }
        }

        //細胞に線のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->li.push_back(i);
        }

        //細胞に点のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->vi.push_back(lp->vi[0]);
          cp->vi.push_back(lp->vi[1]);
          //重複の削除
          sort(cp->vi.begin(), cp->vi.end());
          cp->vi.erase(unique(cp->vi.begin(), cp->vi.end()), cp->vi.end());
        }


        //線が持つ頂点のインデックス・頂点が持つ線のインデックス・点が持つ細胞のインデックスを更新する。
        if (c_vec0 * tmp_v0 > 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[1]) { //この時点で、共有線iのインデックスはないので、これでOK
              p_g->p_l[*it]->vi[1] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[1]) {
              p_g->p_l[*it]->vi[0] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++) {
            if ((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i) {
              vp[0]->li.erase(jt);
              break;
            }
          }

          for (auto it = begin(cp1->li); it != end(cp1->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[0]) {
              p_g->p_l[*it]->vi[1] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[0]) {
              p_g->p_l[*it]->vi[0] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++) {
            if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i) {
              vp[1]->li.erase(jt);
              break;
            }
          }

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          vp[0]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          vp[1]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[0]) {
              p_g->p_l[*it]->vi[1] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[0]) {
              p_g->p_l[*it]->vi[0] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++) {
            if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i) {
              vp[1]->li.erase(jt);
              break;
            }
          }

          for (auto it = begin(cp1->li); it != end(cp1->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[1]) {
              p_g->p_l[*it]->vi[1] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[1]) {
              p_g->p_l[*it]->vi[0] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++) {
            if ((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i) {
              vp[0]->li.erase(jt);
              break;
            }
          }

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          vp[0]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          vp[1]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else {
          std::cout << "Error: recconection was failed." << std::endl;
          exit(0);
        }


        //線が持つ細胞の情報を更新する。
        lp->ci.clear();
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          lp->ci.push_back(ci_tmp[j]);
        }

        /*
        	//線が持つ細胞の情報を更新する。
        	std::vector<int> ci_tmp;
        	for(int j = 0; j < 2; j++){
        		for(auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++){
        			bool flag = false;
        			for(auto jt = begin(lp->ci); jt != end(lp->ci); jt++){
        				if((*it) == (*jt)){
        					break;
        				}
        				if(jt == end(lp->ci)-1){
        					flag = true;
        				}
        			}
        			if(flag == true){
        				ci_tmp.push_back(*it);
        			}
        		}
        	}

        	//細胞から線のインデックスを追加・消去する。
        	//消去
        	for(int j = 0; j < lp->ci.size(); j++){
        		Cellula *cp = p_g->p_c[lp->ci[j]];
        		for(auto it = begin(cp->li); it != end(cp->li); it++){
        			if(*it == i){
        				cp->li.erase(it);
        			}
        		}
        	}
        	//追加
        	for(int j = 0; j < ci_tmp.size(); j++){
        		Cellula *cp = p_g->p_c[ci_tmp[j]];
        		cp.push_back(i);
        	}
        	*/

        //std::cout << "pre_flag = 1 end" << std::endl;
      }
      else if ( pre_flag == 2) {
        p_g->reconnection++;
        cell_rec.push_back(lp->ci[0]);
        cell_rec.push_back(lp->ci[1]);
        //std::cout << "pre_flag = 1 start" << std::endl;
        rc_flag = true;
        //座標の移動
        _vec<double> mid_vec = 0.5 * (vp[0]->loc[0] + vp[1]->loc[0]);
        _vec<double> tmp_v0 = _vec<double>(0.0, 0.0, 0.0);
        _vec<double> tmp_v1 = _vec<double>(0.0, 0.0, 0.0);
        tmp_v0.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[0]->loc[0].y);
        tmp_v0.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[0]->loc[0].x - mid_vec.x);
        tmp_v1.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[1]->loc[0].y);
        tmp_v1.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[1]->loc[0].x - mid_vec.x);

        vp[0]->loc[0] = mid_vec + tmp_v0;
        vp[1]->loc[0] = mid_vec + tmp_v1;

        //線が新しく持つ細胞のインデックス
        std::vector<int> ci_tmp;
        for (int j = 0; j < 2; j++) {
          for (auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++) {
            bool flag = false;
            for (auto jt = begin(lp->ci); jt != end(lp->ci); jt++) {
              if ((*it) == (*jt)) {
                break;
              }
              if (jt == end(lp->ci) - 1) {
                flag = true;
              }
            }
            if (flag == true) {
              ci_tmp.push_back(*it);
            }
          }
        }
        std::cout << "pre_flag = " << pre_flag << " cell size = " << ci_tmp.size() << std::endl;

        //細胞の重心との位置関係で移動後の点がどの位置にあるのか判断する。細胞の変形があまりにも激しいとバグるかもしれない。
        Cellula *cp0 = p_g->p_c[lp->ci[0]];
        Cellula *cp1 = p_g->p_c[lp->ci[1]];
        _vec<double> c_vec0 = cp0->center - mid_vec;
        //もし内積が正（重心と同じ側にvp[0]がある）ならばcp0がvp0の移動先
        if (c_vec0 * tmp_v0 > 0) {
          //点0側の細胞から点1の情報を消す。両方の点についてやる。
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[1]) {
              cp0->vi.erase(it);
              break;
            }
          }
          for (auto it = begin(cp1->vi); it != end(cp1->vi); it++) {
            if (*it == lp->vi[0]) {
              cp1->vi.erase(it);
              break;
            }
          }

        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[0]) {
              cp0->vi.erase(it);
              break;
            }
          }
          for (auto it = begin(cp1->vi); it != end(cp1->vi); it++) {
            if (*it == lp->vi[1]) {
              cp1->vi.erase(it);
              break;
            }
          }

        }
        else {
          std::cout << "Error: recconection was failed. c_vec0*tmp_v0 = " << c_vec0 *tmp_v0 << " " << cp0->center.x << std::endl;
          exit(0);
        }

        //細胞から線のインデックスを消去する。
        for (int j = 0; j < (int)lp->ci.size(); j++) {
          Cellula *cp = p_g->p_c[lp->ci[j]];
          for (auto it = begin(cp->li); it != end(cp->li); it++) {
            if (*it == i) {
              cp->li.erase(it);
              break;
            }
          }
        }

        //細胞に線のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->li.push_back(i);
        }

        //細胞に点のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->vi.push_back(lp->vi[0]);
          cp->vi.push_back(lp->vi[1]);
          //重複の削除
          sort(cp->vi.begin(), cp->vi.end());
          cp->vi.erase(unique(cp->vi.begin(), cp->vi.end()), cp->vi.end());
        }


        //線が持つ頂点のインデックス・頂点が持つ線のインデックス・点が持つ細胞のインデックスを更新する。
        if (c_vec0 * tmp_v0 > 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[1]) { //この時点で、共有線iのインデックスはないので、これでOK
              p_g->p_l[*it]->vi[1] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[1]) {
              p_g->p_l[*it]->vi[0] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }

          for (auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++) {
            /*
            int flag_break = 0;
            for( auto cjt:p_g->p_l[*jt]->ci ){//autoで回すと途中で計算が破綻する。切り替えが上手くできていない。原因不明
            	if( cjt == lp->ci[1] ){
            		vp[0]->li.erase(jt);
            		flag_break = 1;
            		break;
            	}

            }
            if( flag_break == 1 ) break;
            */

            if ( p_g->p_l[*jt]->ci.size() == 2 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i) {
                vp[0]->li.erase(jt);
                break;
              }
            }
            else if ( p_g->p_l[*jt]->ci.size() == 1 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[1] ) && *jt != i) {
                vp[0]->li.erase(jt);
                break;
              }
            }
            else {
              std::cout << "NG" << std::endl;
            }



          }

          for (auto it = begin(cp1->li); it != end(cp1->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[0]) {
              p_g->p_l[*it]->vi[1] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[0]) {
              p_g->p_l[*it]->vi[0] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++) {
            /*
            int flag_break = 0;
            for( auto cjt:p_g->p_l[*jt]->ci ){
            	if( cjt == lp->ci[0] ){
            		vp[1]->li.erase(jt);
            		flag_break = 1;
            		break;
            	}

            }
            if( flag_break == 1 ) break;
            */
            if ( p_g->p_l[*jt]->ci.size() == 2 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i) {
                vp[1]->li.erase(jt);
                break;
              }
            }
            else if ( p_g->p_l[*jt]->ci.size() == 1 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] ) && *jt != i) {
                vp[1]->li.erase(jt);
                break;
              }
            }
            else {
              std::cout << "NG" << std::endl;
            }


          }

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          vp[0]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          vp[1]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[0]) {
              p_g->p_l[*it]->vi[1] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[0]) {
              p_g->p_l[*it]->vi[0] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++) {

            /*
            int flag_break = 0;
            for( auto cjt:p_g->p_l[*jt]->ci ){
            	if( cjt == lp->ci[0] ){
            		vp[1]->li.erase(jt);
            		flag_break = 1;
            		break;
            	}

            }
            if( flag_break == 1 ) break;
            */


            int cis = p_g->p_l[*jt]->ci.size();

            if ( cis == 2 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i) {
                vp[1]->li.erase(jt);
                break;
              }
            }
            else if ( cis == 1 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[0]) && *jt != i) {
                vp[1]->li.erase(jt);
                break;
              }
            }
            else {
              std::cout << "NGx2" << std::endl;
            }

          }

          for (auto it = begin(cp1->li); it != end(cp1->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[1]) {
              p_g->p_l[*it]->vi[1] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[1]) {
              p_g->p_l[*it]->vi[0] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }
          for (auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++) {
            /*
            int flag_break = 0;
            for( auto cjt:p_g->p_l[*jt]->ci ){
            	if( cjt == lp->ci[1] ){
            		vp[0]->li.erase(jt);
            		flag_break = 1;
            		break;
            	}

            }
            if( flag_break == 1 ) break;
            */


            int cis = p_g->p_l[*jt]->ci.size();

            if ( cis == 2 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i) {
                vp[0]->li.erase(jt);
                break;
              }
            }
            else if ( cis == 1 ) {
              if ((p_g->p_l[*jt]->ci[0] == lp->ci[1]) && *jt != i) {
                vp[0]->li.erase(jt);
                break;
              }
            }
            else {
              std::cout << "NGx2" << std::endl;
            }


          }

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          vp[0]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          vp[1]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else {
          std::cout << "Error: recconection was failed." << std::endl;
          exit(0);
        }


        //線が持つ細胞の情報を更新する。
        lp->ci.clear();
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          lp->ci.push_back(ci_tmp[j]);
        }

        /*
        	//線が持つ細胞の情報を更新する。
        	std::vector<int> ci_tmp;
        	for(int j = 0; j < 2; j++){
        		for(auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++){
        			bool flag = false;
        			for(auto jt = begin(lp->ci); jt != end(lp->ci); jt++){
        				if((*it) == (*jt)){
        					break;
        				}
        				if(jt == end(lp->ci)-1){
        					flag = true;
        				}
        			}
        			if(flag == true){
        				ci_tmp.push_back(*it);
        			}
        		}
        	}

        	//細胞から線のインデックスを追加・消去する。
        	//消去
        	for(int j = 0; j < lp->ci.size(); j++){
        		Cellula *cp = p_g->p_c[lp->ci[j]];
        		for(auto it = begin(cp->li); it != end(cp->li); it++){
        			if(*it == i){
        				cp->li.erase(it);
        			}
        		}
        	}
        	//追加
        	for(int j = 0; j < ci_tmp.size(); j++){
        		Cellula *cp = p_g->p_c[ci_tmp[j]];
        		cp.push_back(i);
        	}
        	*/

        //std::cout << "pre_flag = 1 end" << std::endl;
      }
      else if ( pre_flag == 3 ) {  //境界型エッジの切り替え
        p_g->reconnection++;
        cell_rec.push_back(lp->ci[0]);
        //std::cout << "pre_flag = 3 start" << std::endl;
        rc_flag = true;
        //座標の移動
        _vec<double> mid_vec = 0.5 * (vp[0]->loc[0] + vp[1]->loc[0]);
        _vec<double> tmp_v0 = _vec<double>(0.0, 0.0, 0.0);
        _vec<double> tmp_v1 = _vec<double>(0.0, 0.0, 0.0);
        tmp_v0.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[0]->loc[0].y);
        tmp_v0.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[0]->loc[0].x - mid_vec.x);
        tmp_v1.x = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (mid_vec.y - vp[1]->loc[0].y);
        tmp_v1.y = 0.5 * (L_THRESHOLD + L_RECONNECTED) / norm * (vp[1]->loc[0].x - mid_vec.x);

        vp[0]->loc[0] = mid_vec + tmp_v0;
        vp[1]->loc[0] = mid_vec + tmp_v1;

        //線が新しく持つ細胞のインデックス
        std::vector<int> ci_tmp;
        for (int j = 0; j < 2; j++) {
          for (auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++) {
            bool flag = false;
            for (auto jt = begin(lp->ci); jt != end(lp->ci); jt++) {
              if ((*it) == (*jt)) {
                break;
              }
              if (jt == end(lp->ci) - 1) {
                flag = true;
              }
            }
            if (flag == true) {
              ci_tmp.push_back(*it);
            }
          }
        }


        std::cout << "pre_flag = " << pre_flag << " cell size = " << ci_tmp.size() << std::endl;

        //細胞の重心との位置関係で移動後の点がどの位置にあるのか判断する。細胞の変形があまりにも激しいとバグるかもしれない。
        Cellula *cp0 = p_g->p_c[lp->ci[0]];
        //Cellula *cp1 = p_g->p_c[lp->ci[1]];
        _vec<double> c_vec0 = cp0->center - mid_vec;
        //もし内積が正（重心と同じ側にvp[0]がある）ならばcp0がvp0の移動先

        if (c_vec0 * tmp_v0 > 0) {
          //点0側の細胞から点1の情報を消す。両方の点についてやる。
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[1]) {
              cp0->vi.erase(it);
              break;
            }
          }
          /*
          for(auto it = begin(cp1->vi); it != end(cp1->vi); it++){
          	if(*it == lp->vi[0]){
          		cp1->vi.erase(it);
          		break;
          	}
          }*/


        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->vi); it != end(cp0->vi); it++) {
            if (*it == lp->vi[0]) {
              cp0->vi.erase(it);
              break;
            }
          }

          /*
          for(auto it = begin(cp1->vi); it != end(cp1->vi); it++){
          	if(*it == lp->vi[1]){
          		cp1->vi.erase(it);
          		break;
          	}
          }
          */

        }
        else {
          std::cout << "Error: recconection was failed." << std::endl;
          exit(0);
        }



        //細胞から線のインデックスを消去する。
        for (int j = 0; j < (int)lp->ci.size(); j++) {
          Cellula *cp = p_g->p_c[lp->ci[j]];
          for (auto it = begin(cp->li); it != end(cp->li); it++) {
            if (*it == i) {
              cp->li.erase(it);
              break;
            }
          }
        }

        //細胞に線のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->li.push_back(i);
        }

        //細胞に点のインデックスを追加する。
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          Cellula *cp = p_g->p_c[ci_tmp[j]];
          cp->vi.push_back(lp->vi[0]);
          cp->vi.push_back(lp->vi[1]);
          //重複の削除
          sort(cp->vi.begin(), cp->vi.end());
          cp->vi.erase(unique(cp->vi.begin(), cp->vi.end()), cp->vi.end());
        }

        std::cout << "debug: bnd edge reconnect: " << c_vec0 *tmp_v0 << std::endl;

        //線が持つ頂点のインデックス・頂点が持つ線のインデックス・点が持つ細胞のインデックスを更新する。
        if (c_vec0 * tmp_v0 > 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[1]) {
              p_g->p_l[*it]->vi[1] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[1]) {
              p_g->p_l[*it]->vi[0] = lp->vi[0];
              vp[0]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }


          }
          /*
          for(auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++){
          	if((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i){
          		vp[0]->li.erase(jt);
          		break;
          	}
          }
          */
          /*
          for( auto it = begin(vp[0]->li); it != end(vp[0]->li); ++it ){
          if( p_g->p_l[*it]->ci.size() == 1 && *it != i){


          if(p_g->p_l[*it]->vi[1] == lp->vi[0]){
          	p_g->p_l[*it]->vi[1] = lp->vi[1];
          	vp[1]->li.push_back(*it);
          	break;
          }else if(p_g->p_l[*it]->vi[0] == lp->vi[0]){
          	p_g->p_l[*it]->vi[0] = lp->vi[1];
          	vp[1]->li.push_back(*it);
          	break;
          }else{
          	std::cout << "missing 1" << std::endl; exit(0);
          }

          }
          }*/

          int l_v021;
          for ( auto it = begin( vp[0]->li ); it != end( vp[0]->li ); ++it ) {
            if ( p_g->p_l[*it]->ci.size() == 1  && *it != i ) {
              l_v021 = *it;

              if (p_g->p_l[l_v021]->vi[1] == lp->vi[0]) {
                p_g->p_l[l_v021]->vi[1] = lp->vi[1];//ad hocにはここをコメントアウトする
                vp[1]->li.push_back(l_v021);
              }
              else if (p_g->p_l[l_v021]->vi[0] == lp->vi[0]) {
                p_g->p_l[l_v021]->vi[0] = lp->vi[1];//ad hocにはここをコメントアウトする
                vp[1]->li.push_back(l_v021);
              }

              vp[0]->li.erase(it);
              break;
            }
          }

          /*
          if(p_g->p_l[l_v021]->vi[1] == lp->vi[0]){
          			p_g->p_l[l_v021]->vi[1] = lp->vi[1];
          			vp[1]->li.push_back(l_v021);
          }else if(p_g->p_l[l_v021]->vi[0] == lp->vi[0]){
          			p_g->p_l[l_v021]->vi[0] = lp->vi[1];
          			vp[1]->li.push_back(l_v021);
          }
          */

          /* debug */

          /*
          for(auto it = begin(cp1->li); it != end(cp1->li); it++){
          	if(p_g->p_l[*it]->vi[1] == lp->vi[0]){
          		p_g->p_l[*it]->vi[1] = lp->vi[1];
          		vp[1]->li.push_back(*it);
          		break;
          	} else if(p_g->p_l[*it]->vi[0] == lp->vi[0]){
          		p_g->p_l[*it]->vi[0] = lp->vi[1];
          		vp[1]->li.push_back(*it);
          		break;
          	} else{
          		continue;
          	}
          }
          */

          for (auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++) {
            if ((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i) {
              vp[1]->li.erase(jt);
              break;
            }
          }

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          vp[0]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          //vp[1]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else if (c_vec0 * tmp_v0 < 0) {
          for (auto it = begin(cp0->li); it != end(cp0->li); it++) {
            if (p_g->p_l[*it]->vi[1] == lp->vi[0]) {
              p_g->p_l[*it]->vi[1] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else if (p_g->p_l[*it]->vi[0] == lp->vi[0]) {
              p_g->p_l[*it]->vi[0] = lp->vi[1];
              vp[1]->li.push_back(*it);
              break;
            }
            else {
              continue;
            }
          }

          int l_v120;
          for ( auto it = begin( vp[1]->li ); it != end( vp[1]->li ); ++it ) {
            if ( p_g->p_l[*it]->ci.size() == 1  && *it != i ) {
              l_v120 = *it;

              if (p_g->p_l[l_v120]->vi[1] == lp->vi[1]) {
                p_g->p_l[l_v120]->vi[1] = lp->vi[0];
                vp[0]->li.push_back(l_v120);
              }
              else if (p_g->p_l[l_v120]->vi[0] == lp->vi[1]) {
                p_g->p_l[l_v120]->vi[0] = lp->vi[0];
                vp[0]->li.push_back(l_v120);
              }
              vp[1]->li.erase(it);
              break;
            }
          }
          /*
          if(p_g->p_l[l_v120]->vi[1] == lp->vi[1]){
          			p_g->p_l[l_v120]->vi[1] = lp->vi[0];
          			vp[0]->li.push_back(l_v120);
          }else if(p_g->p_l[l_v120]->vi[0] == lp->vi[1]){
          			p_g->p_l[l_v120]->vi[0] = lp->vi[0];
          			vp[0]->li.push_back(l_v120);
          }*/



          /* debug
          for( auto it = begin(vp[1]->li); it != end(vp[1]->li); ++it ){
          	if( p_g->p_l[*it]->ci.size() == 1 && *it != i ){
          		if(p_g->p_l[*it]->vi[1] == lp->vi[1]){
          			p_g->p_l[*it]->vi[1] = lp->vi[0];
          			vp[1]->li.push_back(*it);
          		}else if(p_g->p_l[*it]->vi[0] == lp->vi[1]){
          			p_g->p_l[*it]->vi[0] = lp->vi[0];
          			vp[1]->li.push_back(*it);
          		}

          		vp[1]->li.erase(it);
          		break;
          	}
          }
          */
          /*
          for( auto it = begin( vp[0]->li ); it != end( vp[0]->li ); ++it ){
          	if( p_g->p_l[*it]->ci.size() == 2 ){
          		vp[0]->li.erase(it);
          	}
          }
          */

          /*
          for(auto jt = begin(vp[1]->li); jt != end(vp[1]->li); jt++){//怪しい
          	if((p_g->p_l[*jt]->ci[0] == lp->ci[0] || p_g->p_l[*jt]->ci[1] == lp->ci[0]) && *jt != i){
          		vp[1]->li.erase(jt);
          		break;
          	}
          }
          */
          /*
          for(auto it = begin(cp1->li); it != end(cp1->li); it++){
          	if(p_g->p_l[*it]->vi[1] == lp->vi[1]){
          		p_g->p_l[*it]->vi[1] = lp->vi[0];
          		vp[0]->li.push_back(*it);
          		break;
          	} else if(p_g->p_l[*it]->vi[0] == lp->vi[1]){
          		p_g->p_l[*it]->vi[0] = lp->vi[0];
          		vp[0]->li.push_back(*it);
          		break;
          	} else{
          		continue;
          	}
          }
          */
          /*
          for(auto jt = begin(vp[0]->li); jt != end(vp[0]->li); jt++){
          	if((p_g->p_l[*jt]->ci[0] == lp->ci[1] || p_g->p_l[*jt]->ci[1] == lp->ci[1]) && *jt != i){
          		vp[0]->li.erase(jt);
          		break;
          	}
          }
          */

          //点が持つ細胞のインデックスを更新する。
          vp[0]->ci.clear();
          //vp[0]->ci.push_back(lp->ci[1]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[0]->ci.push_back(ci_tmp[j]);
          }
          vp[1]->ci.clear();
          vp[1]->ci.push_back(lp->ci[0]);
          for (int j = 0; j < (int)ci_tmp.size(); j++) {
            vp[1]->ci.push_back(ci_tmp[j]);
          }

        }
        else {
          std::cout << "Error: recconection was failed." << std::endl;
          exit(0);
        }


        //線が持つ細胞の情報を更新する。
        lp->ci.clear();
        for (int j = 0; j < (int)ci_tmp.size(); j++) {
          lp->ci.push_back(ci_tmp[j]);
        }

        /*
        	//線が持つ細胞の情報を更新する。
        	std::vector<int> ci_tmp;
        	for(int j = 0; j < 2; j++){
        		for(auto it = begin(vp[j]->ci); it != end(vp[j]->ci); it++){
        			bool flag = false;
        			for(auto jt = begin(lp->ci); jt != end(lp->ci); jt++){
        				if((*it) == (*jt)){
        					break;
        				}
        				if(jt == end(lp->ci)-1){
        					flag = true;
        				}
        			}
        			if(flag == true){
        				ci_tmp.push_back(*it);
        			}
        		}
        	}

        	//細胞から線のインデックスを追加・消去する。
        	//消去
        	for(int j = 0; j < lp->ci.size(); j++){
        		Cellula *cp = p_g->p_c[lp->ci[j]];
        		for(auto it = begin(cp->li); it != end(cp->li); it++){
        			if(*it == i){
        				cp->li.erase(it);
        			}
        		}
        	}
        	//追加
        	for(int j = 0; j < ci_tmp.size(); j++){
        		Cellula *cp = p_g->p_c[ci_tmp[j]];
        		cp.push_back(i);
        	}
        	*/
        //std::cout << "pre_flag = 3 end" << std::endl;
      }

      for ( auto cr_i : cell_rec ) {
        p_g->p_c[cr_i]->flag_rec = 1;
      }

    }
  }
  //for debug
  //if(rc_flag == true){
  //	output::outputVTK(&p_g, p_g->step);
  //	exit(0);
  //}

  //1度でも繋ぎ替えが行われたら頂点を各細胞について反時計回りにソートし直す。
  if (rc_flag == true) {
    //sortAntiClockwise(p_g);
    sortCounterClockwise(p_g);
  }

};

} // namespace restructure
