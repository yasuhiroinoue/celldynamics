/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#include <iostream>
#include <vector>
#include <omp.h>

#include "vec.h"
#include "_class_and_variables.h"
//#include <random>

//
_RNG::_RNG(){
	left = 1;
	initf = 0;
}
void _RNG::init_genrand(unsigned long s){

    int j;

    state[0]= s & 0xffffffffUL;

    for (j=1; j<N_MT; j++) {

        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 

        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */

        /* In the previous versions, MSBs of the seed affect   */

        /* only MSBs of the array state[].                        */

        /* 2002/01/09 modified by Makoto Matsumoto             */

        state[j] &= 0xffffffffUL;  /* for >32 bit machines */

    }

    left = 1; initf = 1;

}

void _RNG::next_state(void)
{

    unsigned long *p=state;

    int j;



    /* if init_genrand() has not been called, */

    /* a default initial seed is used         */

    if (initf==0) init_genrand(5489UL);



    left = N_MT;

    next = state;

    

    for (j=N_MT-M_MT+1; --j; p++) 

        *p = p[M_MT] ^ TWIST(p[0], p[1]);



    for (j=M_MT; --j; p++) 

        *p = p[M_MT-N_MT] ^ TWIST(p[0], p[1]);



    *p = p[M_MT-N_MT] ^ TWIST(p[0], state[0]);

}

/* generates a random number on (0,1)-real-interval */

double _RNG::rnd(void)
{// (0,1)の一様乱数

    unsigned long y;



    if (--left == 0) next_state();

    y = *next++;



    /* Tempering */

    y ^= (y >> 11);

    y ^= (y << 7) & 0x9d2c5680UL;

    y ^= (y << 15) & 0xefc60000UL;

    y ^= (y >> 18);



    return ((double)y + 0.5) * (1.0/4294967296.0); 

    /* divided by 2^32 */

}

// ------- 乱数生成 並列計算版 -------------
void Init_Random(unsigned long int _seed) {
  std::cout << "Random seed = " << _seed << std::endl;

  int threads = std::max(1, omp_get_max_threads());

  for (int seed = 0; seed < threads; ++seed) {
    _RNG tmp;
    tmp.init_genrand((seed + _seed) );
    rng.push_back(tmp);
  }
}

double RAND() {
  int id = omp_get_thread_num();
  return rng[id].rnd();
}

double gauss(double t, double u) {

  double ran, gatai;

  do ran = RAND();
  while (ran == 0.0);

  gatai = (double)(sqrt(-2.0 * t * log(ran)) * sin(M_PI * 2.0 * RAND())) + u;

  return (gatai);
}
//

std::vector<_RNG> rng;

//
