
#include "config.h"

/*
 *
 * SCORE FUNCTIONS
 *
 */

int calcMaxScore(Config &config, string &read);
int calcMaxScoreZ(Config &config, vector<int> &offsets);
int scoreZ2(Config &config, vector<long> &locs, int centerIndex, vector<int> &offsets, int num_approx_hits, int num_hits);
int scoreY(Config &config, vector<long> &locs, int center_index, vector<int > &offsets);
int scoreLeft(Config &config, vector<long> &locs, int center_index);
int scoreRight(Config &config, vector<long> &locs, int center_index, int num_hits);
int quickScore(Config &config, vector<long> &locs, int center_index, vector<int> &offsets, int num_approx_hits, int num_hits);
int calcMaxQuickScore(Config &config, vector<int> &offsets, vector<int> keys);
int calcAffineScore(Config &config, vector<long> &loc_array, string &read);
int extendScore(Config &config, string &read, vector<int> &offsets, vector<long> &values, int center_index, vector<long> &loc_array,
			int num_hits, int num_approx_hits, string &whole_genome);
int maxImperfectScore(Config &config, string &read);
