
#include "headers.h"

#pragma once

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

struct Config {

	//	PROGRAM CONST
	int THREAD_NUM;
	bool PRECISE;
	int MULTY_PRECISION;
	int SECOND_PHASE_MULTY_PRECISION;
	int THIRD_PHASE_MULTY_PRECISION;
	bool IS_SECOND_PHASE;
	bool IS_THIRD_PHASE;
	int KEYLEN;
	int KEYLEN2;
	int KEYLEN3;
	int KEYSPACE;
	vector<int> POSITIONS_INDEX;
	int SPLIT_COUNT;
	int FIRST_PHASE_SPLIT_COUNT;
	int SECOND_PHASE_SPLIT_COUNT;
	int TOTAL_BASE_NUM;
	int ALIGNED_BASE_NUM;
	int BUILD_NUMBER;
	int COV_THRES;
	int COV_PADDING;
	int COV_GAPLEN;
	int SECOND_PHASE_MODE;
	bool SECOND_PHASE_PRECISE;
	int SECOND_PHASE_MAX_INDEL2;
	double PRECISION_CUTOFF;
	int TOTAL_BASE_NUM_PREVIOUS;

	int TOTAL_ALIGNED_BASES;

	//	INDEX CONST
	string OUTDIR;
	string INDEX_LOCATION;
	bool OVERWRITE_INDEX;
	long LENGTH_OF_SITES;
	long LENGTH_OF_SIZES;

	//	KEYS CONST
	int MAX_DESIRED_KEYS;
	int MIN_KEYS_DESIRED;
	float MIN_KEY_DENSITY;
	float KEY_DENSITY;
	float MAX_KEY_DENSITY;
	//	UPPER LIMIT FOR KEY OCCURANCE
	int MAX_LEN;

	//	HITS CUTOFF CONST
	int MIN_APPROX_HITS_TO_KEEP;
	int HIT_REDUCTION_DIV;
	int MAX_HITS_REDUCTION2;
	int MAXIMUM_MAX_HITS_REDUCTION;

	//	ALIGNER CONST
	int Y_SCORE_MULT;
	int Z_SCORE_MULT;
	int BASE_SCORE;
	float MIN_QSCORE_MULT;
	float MIN_QSCORE_MULT2;
	float MIN_SCORE_MULT;
	float DYNAMIC_QSCORE_THRESH_PERFECT;
	float DYNAMIC_QSCORE_THRESH;
	float DYNAMIC_SCORE_THRESH;
	int MAX_INDEL;
	int MAX_INDEL2;
	int THIRD_PHASE_MAX_INDEL2;
	float PRESCAN_QSCORE_THRESH;

	int SLOW_ALIGN_PADDING;

	//	FILL UNLIMITED CONST
	long TIMEMASK;
	long SCOREMASK;
	long MAX_TIME;
	long BARRIER_I1;
	long BARRIER_D1;
	long POINTSoff_MATCH2;
	long POINTSoff_MATCH;
	long POINTSoff_SUBR;
	long POINTSoff_SUB;
	long POINTSoff_NOCALL;
	long POINTSoff_DEL;
	long POINTSoff_DEL2;
	long POINTSoff_DEL3;
	long POINTSoff_DEL4;
	long POINTSoff_DEL5;
	long POINTSoff_DEL_REF_N;
	long POINTSoff_GAP;
	long POINTSoff_INS;
	long BADoff;
	int SCOREOFFSET;

	//	SCORE CONST
	int SCOREZ_1KEY;
	int INDEL_PENALTY;
	int INDEL_PENALTY_MULT;
	int MAX_PENALTY_FOR_MISALIGNED_HIT;

	//	MATCH STRING CONST
	int MODE_DEL;
	int MODE_INS;
	int MODE_MS;
	int MODE_SUB;
	int POINTS_MATCH;
	int POINTS_MATCH2;
	int POINTS_INS;
	int POINTS_SUB;
	int POINTS_SUB2;
	int POINTS_SUB3;
	int POINTS_NOREF;
	int MASK5;
	int TIMESLIP;
	int LIMIT_FOR_COST_5;
	int LIMIT_FOR_COST_4;
	int LIMIT_FOR_COST_3;
	int POINTS_NOCALL;
	int POINTS_DEL;
	int POINTS_DEL2;
	int POINTS_DEL3;
	int POINTS_DEL4;
	int POINTS_DEL5;

	//	GAP CONST
	int MINGAP;
	int GAPLEN;
	int GAPBUFFER2;
	int GAPBUFFER;
	char GAPC;
	int POINTS_GAP;

	int MAX_SUBSUMPTION_LENGTH;

	Config() {
		//	PROGRAM CONST
		THREAD_NUM = 4;
		PRECISE = true;
		MULTY_PRECISION = 0;
		SECOND_PHASE_MULTY_PRECISION = 0;
		THIRD_PHASE_MULTY_PRECISION = 0;
		IS_SECOND_PHASE = false;
		IS_THIRD_PHASE = false;
		KEYLEN = 13;
		KEYLEN2 = 13;
		KEYLEN3 = 13;
		KEYSPACE = pow(2, 2*KEYLEN);
		SPLIT_COUNT = 1;
		FIRST_PHASE_SPLIT_COUNT = SPLIT_COUNT;
		SECOND_PHASE_SPLIT_COUNT = 1;
		TOTAL_BASE_NUM = 0;
		ALIGNED_BASE_NUM = 0;
		BUILD_NUMBER = 1;
		COV_THRES = 1;
		COV_PADDING = 100;
		COV_GAPLEN = 32000;
		SECOND_PHASE_MODE = 0;
		SECOND_PHASE_PRECISE  = true;
		SECOND_PHASE_MAX_INDEL2 = 32000;
		PRECISION_CUTOFF = 0.9;
		TOTAL_BASE_NUM_PREVIOUS = TOTAL_BASE_NUM;

		TOTAL_ALIGNED_BASES = 0;

		//	INDEX CONST
		OVERWRITE_INDEX = false;
		LENGTH_OF_SITES = 0;
		LENGTH_OF_SIZES = 0;

		//	KEYS CONST
		MAX_DESIRED_KEYS = 15;
		MIN_KEYS_DESIRED = 2;
		MIN_KEY_DENSITY = 1.5;
		KEY_DENSITY = 1.9;
		MAX_KEY_DENSITY = 3.0;
		//	UPPER LIMIT FOR KEY OCCURANCE
		MAX_LEN = 10998;

		//	HITS CUTOFF CONST
		MIN_APPROX_HITS_TO_KEEP = 1;
		HIT_REDUCTION_DIV = 4;
		MAX_HITS_REDUCTION2 = 3;
		MAXIMUM_MAX_HITS_REDUCTION = 5;

		//	ALIGNER CONST
		Y_SCORE_MULT = 10;
		Z_SCORE_MULT = 25;
		BASE_SCORE = 100*KEYLEN;
		MIN_QSCORE_MULT = 0.005;
		MIN_QSCORE_MULT2 = 0.005;
		MIN_SCORE_MULT = 0.02;
		DYNAMIC_QSCORE_THRESH_PERFECT = 0.8;
		DYNAMIC_QSCORE_THRESH = 0.6;
		DYNAMIC_SCORE_THRESH = 0.64;
		MAX_INDEL = 16000;
		MAX_INDEL2 = 32000;
		THIRD_PHASE_MAX_INDEL2 = MAX_INDEL2;
		PRESCAN_QSCORE_THRESH = 0.6 * 0.95;

		SLOW_ALIGN_PADDING = 6;

		//	FILL UNLIMITED CONST
		TIMEMASK = 2047;
		SCOREMASK = -2048;
		MAX_TIME = 2047;
		BARRIER_I1 = 2;
		BARRIER_D1 = 3;
		POINTSoff_MATCH2 = 204800;
		POINTSoff_MATCH = 143360;
		POINTSoff_SUBR = -301056;
		POINTSoff_SUB = -260096;
		POINTSoff_NOCALL = 0;
		POINTSoff_DEL = -966656;
		POINTSoff_DEL2 = -67584;
		POINTSoff_DEL3 = -18432;
		POINTSoff_DEL4 = -2048;
		POINTSoff_DEL5 = -2048;
		POINTSoff_DEL_REF_N = -20480;
		POINTSoff_GAP = -4096;
		POINTSoff_INS = -808960;
		BADoff = -2143387648;
		SCOREOFFSET = 11;

		//	SCORE CONST
		SCOREZ_1KEY = 260;
		INDEL_PENALTY = 649;
		INDEL_PENALTY_MULT = 20;
		MAX_PENALTY_FOR_MISALIGNED_HIT = 1137;

		//	MATCH STRING CONST
		MODE_DEL = 1;
		MODE_INS = 2;
		MODE_MS = 0;
		MODE_SUB = 3;
		POINTS_MATCH = 70;
		POINTS_MATCH2 = 100;
		POINTS_INS = -395;
		POINTS_SUB = -127;
		POINTS_SUB2 = -51;
		POINTS_SUB3 = -25;
		POINTS_NOREF = 0;
		MASK5 = 3;
		TIMESLIP = 4;
		LIMIT_FOR_COST_5 = 80;
		LIMIT_FOR_COST_4 = 20;
		LIMIT_FOR_COST_3 = 5;
		POINTS_NOCALL = 0;
		POINTS_DEL = -472;
		POINTS_DEL2 = -33;
		POINTS_DEL3 = -9;
		POINTS_DEL4 = -1;
		POINTS_DEL5 = -1;

		//	GAP CONST
		MINGAP = 256;
		GAPLEN = 128;
		GAPBUFFER2 = 128;
		GAPBUFFER = 64;
		GAPC = '-';
		POINTS_GAP = -2;

		MAX_SUBSUMPTION_LENGTH = MAX_INDEL2;
	}
};

struct Result {
	int br;
	long start;
	long stop;
	int score;
	int maxScore;
	string matchString;
	vector<long> gapArray;
	Result(int br_, long start_, long stop_, int score_, int maxScore_) {
		br = br_;
		start = start_;
		stop = stop_;
		score = score_;
		maxScore = maxScore_;
	}
};

struct Read {
	int br;
	string content;
	Read(int br_, string content_) {
		br = br_;
		content = content_;
	}
};

struct FastaRead {
	string name;
	string read;
	FastaRead(string name_, string read_) {
		name = name_;
		read = read_;
	}
};

struct ReferenceSegment {
	int br;
	string name;
	int start;
	ReferenceSegment(int br_, string name_, int start_) {
		br = br_;
		name = name_;
		start = start_;
	}
};

struct ThreadData3 {
	int thread_id;
	long *sizes;
	long *sites;
	vector<Read> *reads;
	int start;
	int stop;
	vector<vector<Result> > *databaseResults;
	Config *config;
	string *whole_genome;
	vector<Result> *results;
	vector<Read> *aligned_reads;
	vector<Read> *unaligned_reads;
	ThreadData3() {
	}
	ThreadData3(int thread_id_, long *sizes_, long *sites_, vector<Read> *reads_, int start_, int stop_, Config *config_, string *whole_genome_, vector<Result> *results_, vector<Read> *aligned_reads_, vector<Read> *unaligned_reads_, vector<vector<Result> > *databaseResults_) {
		thread_id = thread_id_;
		sizes = sizes_;
		sites = sites_;
		reads = reads_;
		start = start_;
		stop = stop_;
		config = config_;
		databaseResults = databaseResults_;
		whole_genome = whole_genome_;
		results = results_;
		aligned_reads = aligned_reads_;
		unaligned_reads = unaligned_reads_;
	}
};

struct SiteScore {
	long start;
	long stop;
	int score;
	int hits;
	bool perfect;
	int strand;
	vector<long> gapArray;
	SiteScore() {
		start = -1;
		stop = -1;
		score = -1;
		hits = -1;
		perfect = false;
		strand = -1;
	}
	SiteScore(long start_, long stop_, int score_, int hits_, bool perfect_, int strand_, vector<long> gapArray_) {
		start = start_;
		stop = stop_;
		score = score_;
		hits = hits_;
		perfect = perfect_;
		strand = strand_;
		gapArray = gapArray_;
	}
};
