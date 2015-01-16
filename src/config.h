#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h>

using namespace std;

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

static string outdir;
static string index_location;
static bool overwrite_index = false;
static int thread_num = 4;

static int length_of_sites;
static int length_of_sizes;

static bool precise = true;

static int maxDesiredKeys = 15;
static int minKeysDesired = 2;
static float minKeyDensity = 1.5;
static float keyDensity = 1.9;
static float maxKeyDensity = 3.0;

static int POINTS_INS = -395;

static int MAX_LEN = 10998;

static int MIN_APPROX_HITS_TO_KEEP = 1;

static int HIT_REDUCTION_DIV = 4;
static int MAX_HITS_REDUCTION2 = 3;
//static int MAXIMUM_MAX_HITS_REDUCTION = 5;
static int MAXIMUM_MAX_HITS_REDUCTION = 6;

static bool drugi = false;

static vector<int> positions_index;

static int split_count = 1;
static int total_base_num = 0;
static int aligned_base_num = 0;

static int KEYLEN = 13;
static int keyspace = pow(2, 2*KEYLEN);

static int build_number = 1;

static float MIN_QSCORE_MULT2 = 0.005;
//static float MIN_QSCORE_MULT2 = 0.1;
//static float PRESCAN_QSCORE_THRESH = 0.57;
static float PRESCAN_QSCORE_THRESH = 0.6 * 0.95;

static int MODE_INS = 2;
static int MODE_MS = 0;
static int MODE_DEL = 1;
static int MODE_SUB = 3;
static int POINTS_SUB2 = -51;
static int POINTS_SUB3 = -25;
static int POINTS_NOREF = 0;
static int GAPBUFFER = 64;
static char GAPC = '-';

static int SLOW_ALIGN_PADDING = 4;
static int MARIC_PADDING = 6;

static long POINTSoff_MATCH2 = 204800;
static long POINTSoff_MATCH = 143360;

static long BARRIER_I1 = 2;
static long BARRIER_D1 = 3;

static long TIMEMASK = 2047;
static long SCOREMASK = -2048;
static long MAX_TIME = 2047;

static long POINTSoff_SUBR = -301056;
static long POINTSoff_SUB = -260096;
static long POINTSoff_NOCALL = 0;

static long POINTSoff_DEL = -966656;
static long POINTSoff_DEL2 = -67584;
static long POINTSoff_DEL3 = -18432;
static long POINTSoff_DEL4 = -2048;
static long POINTSoff_DEL5 = -2048;
static long POINTSoff_DEL_REF_N = -20480;
static long POINTSoff_GAP = -4096;

static long POINTSoff_INS = -808960;
static int SCOREOFFSET = 11;

static int POINTS_MATCH = 70;
static int POINTS_MATCH2 = 100;
//static int Z_SCORE_MULT = 20;
static int Z_SCORE_MULT = 25;
static int BASE_SCORE = 1300;
static int Y_SCORE_MULT = 10;
static float MIN_QSCORE_MULT = 0.005;
//static float MIN_QSCORE_MULT = 0.025;
static float MIN_SCORE_MULT = 0.02;
//static float MIN_SCORE_MULT = 0.15;
static float DYNAMIC_QSCORE_THRESH_PERFECT = 0.8;
//static float DYNAMIC_QSCORE_THRESH_PERFECT = 0.8;
//static float DYNAMIC_QSCORE_THRESH = 0.6;
static float DYNAMIC_QSCORE_THRESH = 0.6;
static float DYNAMIC_SCORE_THRESH = 0.64;
//static float DYNAMIC_SCORE_THRESH = 0.84;
static int MAX_INDEL = 16000;
static int MAX_INDEL2 = 16000;

// score constants
static int SCOREZ_1KEY = 260;
static int INDEL_PENALTY = 649;
static int INDEL_PENALTY_MULT = 20;
static int MAX_PENALTY_FOR_MISALIGNED_HIT = 1137;
static int POINTS_DEL = -472;
static int POINTS_NOCALL = 0;
static int MASK5 = 3;
static int TIMESLIP = 4;
static int LIMIT_FOR_COST_5 = 80;
static int LIMIT_FOR_COST_4 = 20;
static int POINTS_DEL4 = -1;
static int LIMIT_FOR_COST_3 = 5;
static int POINTS_DEL3 = -9;
static int POINTS_DEL2 = -33;
static int POINTS_SUB = -127;
static int POINTS_DEL5 = -1;
static int MINGAP = 256;
static int GAPLEN = 128;
static int GAPBUFFER2 = 128;
static int POINTS_GAP = -2;

static long BADoff = -2143387648;

static int MAX_SUBSUMPTION_LENGTH = 32000;

