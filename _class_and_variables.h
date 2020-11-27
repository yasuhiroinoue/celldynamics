/*****************/
// 2D Vertex Model
// Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp)
/*****************/
#ifndef _CLASS_VAR_H
#define _CLASS_VAR_H

#include <vector>
#include <omp.h>

#include "vec.h"
//#include <random>

//メルセンヌ・ツイスタ乱数用 ここから //
#define N_MT 624
#define M_MT 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))


class _RNG{
private:
/* Period parameters */  

	unsigned long state[N_MT]; /* the array for the state vector  */
	int left;
	int initf;
	unsigned long *next;
	
public:
	/* initializes state[N] with a seed */
	_RNG();
	void next_state(void);
	void init_genrand(unsigned long s);
	double rnd(void);
};

void Init_Random(unsigned long ing);
double RAND(void);
double gauss(double, double);

extern std::vector<_RNG> rng;//乱数

//メルセンヌ・ツイスタ乱数用　ここまで //


/*
class RNG
{//掲示板からパクってきた http://stackoverflow.com/questions/15918758/how-to-make-each-thread-use-its-own-rng-in-c11
public:
    typedef std::mt19937 Engine;
    typedef std::uniform_real_distribution<double> Distribution;

    RNG() : engines(), distribution(0.0, 1.0)
    {
	//入力seedを利用するように変更
	unsigned long int _seed = 20110412;
	char fna[100] = "seed.txt";
	FILE* fptr = fopen(fna,"r");
	fscanf(fptr,"%d",&_seed);
	cout << "Random seed = " << _seed << endl;
	fclose(fptr);
	
        int threads = std::max(1, omp_get_max_threads());
        for(int seed = 0; seed < threads; ++seed)
        {
            engines.push_back(Engine((seed+_seed)));
        }
    }

    double operator()()
    {
        int id = omp_get_thread_num();
        return distribution(engines[id]);
    }

    std::vector<Engine> engines;
    Distribution distribution;
};

RNG rng;

#define RAND rng()
*/


#endif

