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

struct Statistic {
	int total_reads;
	int aligned_reads;
	int unaligned_reads;
	int start15;
	int stop15;
	int start_and_stop15;
	int start_or_stop15;
	int start30;
	int stop30;
	int start_and_stop30;
	int start_or_stop30;
	int start50;
	int stop50;
	int start_and_stop50;
	int start_or_stop50;
	int start100;
	int stop100;
	int start_and_stop100;
	int start_or_stop100;
	int total_exons;
	int overlaping_exons;
	int found_exons;
	int total_bases;
	int aligned_bases;
	int covered_bases;
	int finded_bases;
	Statistic() {
	}
	Statistic(int num_reads) {
		total_reads = num_reads;
		aligned_reads = 0;
		unaligned_reads = 0;
		start15 = 0;
		stop15 = 0;
		start_and_stop15 = 0;
		start_or_stop15 = 0;
		start30 = 0;
		stop30 = 0;
		start_and_stop30 = 0;
		start_or_stop30 = 0;
		start50 = 0;
		stop50 = 0;
		start_and_stop50 = 0;
		start_or_stop50 = 0;
		start100 = 0;
		stop100 = 0;
		start_and_stop100 = 0;
		start_or_stop100 = 0;
		total_exons = 0;
		overlaping_exons = 0;
		found_exons = 0;
		total_bases = 0;
		aligned_bases = 0;
		covered_bases = 0;
		finded_bases = 0;
	}
};


struct Info {
	int read_number;
	int aligned_read_number;
	int genome_size;
	double total_execution_time;
	bool first_phase_precision;
	bool second_phase_precision;
	bool third_phase_precision;
	int first_phase_readlen;
	int second_phase_readlen;
	int third_phase_readlen;
	int first_phase_keylen;
	int second_phase_keylen;
	int third_phase_keylen;
	int first_phase_max_indel;
	int second_phase_max_indel;
	int third_phase_max_indel;
	int first_phase_multy_precision;
	int second_phase_multy_precision;
	int third_phase_multy_precision;
	int coverage_threshold;
	int coverage_gaplen;
	int coverage_padding;
	int second_phase_mode;
	int new_ref_size;
	int number_of_threads;
	int better;
	int worse;
	int same;
	int unfound;
	double index_generation_time;
	double first_phase_time;
	double coverage_time;
	double new_ref_time;
	double second_ref_indexing_time;
	double second_align_time;
	double third_ref_indexing_time;
	double unprecised_reads_time;
	int memory_after_first_align;
	int memory_after_second_align;
	int memory_after_third_align;
	int memory_at_the_end;
	double first_align_time2;
	double second_align_time2;
	double third_align_time2;
	Info() {
		read_number = 0;
		aligned_read_number = 0;
		genome_size = 0;
		total_execution_time = 0;
		first_phase_precision = true;
		second_phase_precision = true;
		third_phase_precision = true;
		first_phase_readlen = 0;
		second_phase_readlen = 0;
		third_phase_readlen = 0;
		first_phase_keylen = 13;
		second_phase_keylen = 13;
		third_phase_keylen = 13;
		first_phase_max_indel = 32000;
		second_phase_max_indel = 32000;
		third_phase_max_indel = 32000;
		first_phase_multy_precision = 0;
		second_phase_multy_precision = 0;
		third_phase_multy_precision = 0;
		coverage_threshold = 0;
		coverage_gaplen = 32000;
		coverage_padding = 100;
		second_phase_mode = 0;
		new_ref_size = 0;
		number_of_threads = 4;
		index_generation_time = 0.0;
		first_phase_time = 0.0;
		coverage_time = 0.0;
		new_ref_time = 0.0;
		second_ref_indexing_time = 0.0;
		third_phase_multy_precision = 0;
		second_align_time = 0.0;
		unprecised_reads_time = 0.0;
		third_ref_indexing_time = 0.0;
		better = 0;
		worse = 0;
		same = 0;
		unfound = 0;
		memory_after_first_align = 0;
		memory_after_second_align = 0;
		memory_after_third_align = 0;
		memory_at_the_end = 0;
		first_align_time2  = 0;
		second_align_time2 = 0;
		third_align_time2 = 0;
	}
};

static Info info;
static Statistic global_stat;
static map<int, int> ssstat;
static map<int, int> ststat;

//	PROGRAM CONST
static int THREAD_NUM = 4;
static bool PRECISE = true;
static int MULTY_PRECISION = 0;
static int SECOND_PHASE_MULTY_PRECISION = 0;
static int THIRD_PHASE_MULTY_PRECISION = 0;
static bool IS_SECOND_PHASE = false;
static bool IS_THIRD_PHASE = false;
static int KEYLEN = 13;
static int KEYLEN2 = 13;
static int KEYLEN3 = 13;
static int KEYSPACE = pow(2, 2*KEYLEN);
static vector<int> POSITIONS_INDEX;
static int SPLIT_COUNT = 1;
static int FIRST_PHASE_SPLIT_COUNT = SPLIT_COUNT;
static int SECOND_PHASE_SPLIT_COUNT = 1;
static int TOTAL_BASE_NUM = 0;
static int ALIGNED_BASE_NUM = 0;
static int BUILD_NUMBER = 1;
static int COV_THRES = 1;
static int COV_PADDING = 100;
static int COV_GAPLEN = 32000;
static int SECOND_PHASE_MODE = 0;
static bool SECOND_PHASE_PRECISE  = true;
static int SECOND_PHASE_MAX_INDEL2 = 32000;
static double PRECISION_CUTOFF = 0.9;
static int TOTAL_BASE_NUM_PREVIOUS = TOTAL_BASE_NUM;

//	INDEX CONST
static string OUTDIR;
static string INDEX_LOCATION;
static bool OVERWRITE_INDEX = false;
static int LENGTH_OF_SITES;
static int LENGTH_OF_SIZES;

//	KEYS CONST
static int MAX_DESIRED_KEYS = 15;
static int MIN_KEYS_DESIRED = 2;
static float MIN_KEY_DENSITY = 1.5;
static float KEY_DENSITY = 1.9;
static float MAX_KEY_DENSITY = 3.0;
//	UPPER LIMIT FOR KEY OCCURANCE
static int MAX_LEN = 10998;

//	HITS CUTOFF CONST
static int MIN_APPROX_HITS_TO_KEEP = 1;
static int HIT_REDUCTION_DIV = 4;
static int MAX_HITS_REDUCTION2 = 3;
static int MAXIMUM_MAX_HITS_REDUCTION = 5;

//	ALIGNER CONST
static int Y_SCORE_MULT = 10;
static int Z_SCORE_MULT = 25;
static int BASE_SCORE = 100*KEYLEN;
static float MIN_QSCORE_MULT = 0.005;
static float MIN_QSCORE_MULT2 = 0.005;
static float MIN_SCORE_MULT = 0.02;
static float DYNAMIC_QSCORE_THRESH_PERFECT = 0.8;
static float DYNAMIC_QSCORE_THRESH = 0.6;
static float DYNAMIC_SCORE_THRESH = 0.64;
static int MAX_INDEL = 16000;
static int MAX_INDEL2 = 32000;
static int THIRD_PHASE_MAX_INDEL2 = MAX_INDEL2;
static float PRESCAN_QSCORE_THRESH = 0.6 * 0.95;

static int SLOW_ALIGN_PADDING = 6;

//	FILL UNLIMITED CONST
static long TIMEMASK = 2047;
static long SCOREMASK = -2048;
static long MAX_TIME = 2047;
static long BARRIER_I1 = 2;
static long BARRIER_D1 = 3;
static long POINTSoff_MATCH2 = 204800;
static long POINTSoff_MATCH = 143360;
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
static long BADoff = -2143387648;
static int SCOREOFFSET = 11;

//	SCORE CONST
static int SCOREZ_1KEY = 260;
static int INDEL_PENALTY = 649;
static int INDEL_PENALTY_MULT = 20;
static int MAX_PENALTY_FOR_MISALIGNED_HIT = 1137;

//	MATCH STRING CONST
static int MODE_DEL = 1;
static int MODE_INS = 2;
static int MODE_MS = 0;
static int MODE_SUB = 3;
static int POINTS_MATCH = 70;
static int POINTS_MATCH2 = 100;
static int POINTS_INS = -395;
static int POINTS_SUB = -127;
static int POINTS_SUB2 = -51;
static int POINTS_SUB3 = -25;
static int POINTS_NOREF = 0;
static int MASK5 = 3;
static int TIMESLIP = 4;
static int LIMIT_FOR_COST_5 = 80;
static int LIMIT_FOR_COST_4 = 20;
static int LIMIT_FOR_COST_3 = 5;
static int POINTS_NOCALL = 0;
static int POINTS_DEL = -472;
static int POINTS_DEL2 = -33;
static int POINTS_DEL3 = -9;
static int POINTS_DEL4 = -1;
static int POINTS_DEL5 = -1;
//	GAP CONST
static int MINGAP = 256;
static int GAPLEN = 128;
static int GAPBUFFER2 = 128;
static int GAPBUFFER = 64;
static char GAPC = '-';
static int POINTS_GAP = -2;

static int MAX_SUBSUMPTION_LENGTH = MAX_INDEL2;

