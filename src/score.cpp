

#include "headers.h"
#include "config.h"
#include "const.h"


static string tmp = "testout";
static ofstream out(tmp.c_str());

/*
 *
 * SCORE FUNCTIONS
 *
 */
int calcMaxScore(Config &config, string &read){
	return config.POINTS_MATCH + (read.size()-1) * config.POINTS_MATCH2;
}

int calcMaxScoreZ(Config &config, vector<int> &offsets){
	int score = 0;
	int a0 = -1;
	int b0 = -1;

	for(unsigned int i = 0; i < offsets.size(); i++){
		int a=offsets[i];
		if(b0 < a){
			score += b0-a0;
			a0 = a;
		}
		b0 = a + config.KEYLEN;
	}
	score += b0 - a0;
	return score * config.Z_SCORE_MULT;
}

int scoreZ2(Config &config, vector<long> &locs, int centerIndex, vector<int> &offsets, int num_approx_hits, int num_hits){

	if(num_approx_hits == 1){
		return config.SCOREZ_1KEY;
	}

	long center = locs[centerIndex];
	long maxLoc = center + config.MAX_INDEL2;
	long minLoc = max((long)0, center - config.MAX_INDEL);
	int score=0;
	int a0 = -1;
	int b0 = -1;

	for(int i = 0; i < num_hits; i++){
		long loc = locs[i];
		if(loc >= minLoc && loc <= maxLoc){
			int a = offsets[i];
			if(b0 < a){
				score += b0 - a0;
				a0 = a;
			}
			b0 = a + config.KEYLEN;
		}
	}
	score += b0 - a0;
	score = score * config.Z_SCORE_MULT;
	return score;
}

int scoreY(Config &config, vector<long> &locs, int center_index, vector<int > &offsets) {
		long center = locs[center_index];

		int rightIndex = -1;
		for(int i = offsets.size() - 1; rightIndex < center_index; i--){
			if(locs[i] == center){
				rightIndex = i;
			}
		}
		return offsets[rightIndex] - offsets[center_index];
}

int scoreLeft(Config &config, vector<long> &locs, int center_index) {
	int score=0;
	long prev;
	long loc = locs[center_index];

	for(int i = center_index-1; i >= 0; i--){

		if(locs[i] >= 0){
			prev = loc;
			loc = locs[i];

			int offset = abs(loc - prev);
			if(offset <= config.MAX_INDEL){
				score += config.BASE_SCORE;
				if(offset != 0){
					int penalty = min(config.INDEL_PENALTY+config.INDEL_PENALTY_MULT*offset, config.MAX_PENALTY_FOR_MISALIGNED_HIT);
					score-=penalty;
				}
			}else{
				loc=prev;
			}
		}
	}
	return score;
}

int scoreRight(Config &config, vector<long> &locs, int center_index, int num_hits) {
	int score=0;
	long prev;
	long loc = locs[center_index];

	for(int i = center_index+1; i < num_hits; i++){

		if(locs[i] >= 0){
			prev = loc;
			loc = locs[i];
			int offset = abs(loc - prev);

			if(offset <= config.MAX_INDEL){
				score += config.BASE_SCORE;

				if(offset != 0){
					int penalty = min(config.INDEL_PENALTY + config.INDEL_PENALTY_MULT * offset, config.MAX_PENALTY_FOR_MISALIGNED_HIT);
					score -= penalty;
				}
			}else{
				loc = prev;
			}
		}
	}
	return score;
}

int quickScore(Config &config, vector<long> &locs, int center_index, vector<int> &offsets, int num_approx_hits, int num_hits) {

	if(num_approx_hits == 1) {
		return config.BASE_SCORE;
	}

	int x = config.BASE_SCORE + scoreLeft(config, locs, center_index)+
		scoreRight(config, locs, center_index, num_hits) - center_index;

	int y = config.Y_SCORE_MULT * scoreY(config, locs, center_index, offsets);
	return x+y;
}

int calcMaxQuickScore(Config &config, vector<int> &offsets, vector<int> keys) {
	int x = keys.size() * config.BASE_SCORE;
	int y = config.Y_SCORE_MULT * (offsets[offsets.size()-1] - offsets[0]);
	x += calcMaxScoreZ(config, offsets);
	return x+y;
}

int calcAffineScore(Config &config, vector<long> &loc_array, string &read) {
	int score = 0;
	long lastLoc = -3;
	long lastValue = -1;
	int timeInMode = 0;

	for(unsigned int i = 0; i < loc_array.size(); i++){
		long loc = loc_array[i];

		if(loc > 0) {//match
			if(loc == lastValue){//contiguous match
				score += (config.POINTS_MATCH2);
			}else if(loc == lastLoc || lastLoc < 0) {//match after a sub, or first match
				score += (config.POINTS_MATCH);
			}else if(loc < lastLoc) {//deletion
				score += (config.POINTS_MATCH);
				score += config.POINTS_DEL;
				long dif = lastLoc-loc+1;
				if(dif > config.MINGAP){
					long rem = dif % config.GAPLEN;
					long div = (dif - config.GAPBUFFER2) / config.GAPLEN;
					score += (div * config.POINTS_GAP);
					dif = rem + config.GAPBUFFER2;
				}
				if(dif>config.LIMIT_FOR_COST_5){
					score+=((dif-config.LIMIT_FOR_COST_5+config.MASK5)/config.TIMESLIP)*config.POINTS_DEL5;
					dif=config.LIMIT_FOR_COST_5;
				}
				if(dif>config.LIMIT_FOR_COST_4){
					score+=(dif-config.LIMIT_FOR_COST_4)*config.POINTS_DEL4;
					dif=config.LIMIT_FOR_COST_4;
				}
				if(dif>config.LIMIT_FOR_COST_3){
					score+=(dif-config.LIMIT_FOR_COST_3)*config.POINTS_DEL3;
					dif=config.LIMIT_FOR_COST_3;
				}
				if(dif>1){
					score+=(dif-1)*config.POINTS_DEL2;
				}
				timeInMode=1;
			}else if(loc > lastLoc) {//insertion
				score += (config.POINTS_MATCH + POINTS_INS_ARRAY_C[min(loc-lastLoc, (long)5)]);

				timeInMode=1;
			}
			lastLoc=loc;
		}else if(loc==-1){//substitution
			if(lastValue < 0 && timeInMode > 0){//contiguous
				timeInMode++;
				long temp = POINTS_SUB_ARRAY[timeInMode];
				score += temp;
			}else{
				score += config.POINTS_SUB;
				timeInMode=1;
			}
		}else{
			timeInMode=0;
			score += config.POINTS_NOCALL;
		}
		lastValue=loc;
	}
	//cout << score << "\n";


	int var = 0;
	for(int i = 0; i < 30; i++) {
		if(loc_array[i] < 0) {
			var++;
		}
	}
	if(var > 10) {
		return 0;
	}
	var = 0;
	for(int i = loc_array.size()-1; i > loc_array.size()-30; i--) {
		if(loc_array[i] < 0) {
			var++;
		}
	}
	if(var > 10) {
		return 0;
	}


	return score;
}

int extendScore(Config &config, string &read, vector<int> &offsets, vector<long> &values, int center_index, vector<long> &loc_array, int num_hits, int num_approx_hits, string &whole_genome) {

	long center_loc = values[center_index];
	long min_loc = center_loc - config.MAX_INDEL;
	long max_loc = center_loc + config.MAX_INDEL2;

	loc_array.clear();
	for(unsigned int i = 0; i < read.size(); i++) {
		loc_array.push_back(-1);
	}

	//First fill in reverse
	for(int i = 0, keynum = 0; i < num_hits; i++){
		long value=values[i];
		if(value >= min_loc && value <= max_loc){
			long refbase = value;
			keynum++;
			long callbase = offsets[i];
			int misses=0;

			for(long cloc = callbase+config.KEYLEN-1, rloc = refbase+cloc; cloc >= 0 && rloc >= 0 && rloc < (long) whole_genome.size(); cloc--, rloc--){
				long old = loc_array[cloc];
				if(old == refbase){
					break;
				} //Already filled with present value
				if(misses > 0 && old >= 0){
					break;
				} //Already filled with something that has no errors
				char c = toupper(read[cloc]);
				char r = toupper(whole_genome[rloc]);

				if(c == r){
					if(old < 0 || refbase == center_loc){ //If the cell is empty or this key corresponds to center
						loc_array[cloc]=refbase;
					}
				}else{
					misses++;
					//Only extends first key all the way back.  Others stop at the first error.
					if(old >= 0 || keynum > 1){
						break;
					}
				}
			}
		}
	}

	//Then fill forward
	for(int i = 0; i < num_hits; i++){
		long value = values[i];
		if(value >= min_loc && value <= max_loc){
			long refbase = value;
			long callbase = offsets[i];
			int misses=0;
			for(long cloc = callbase+config.KEYLEN, rloc = refbase+cloc; cloc < (int) read.size() && rloc < (long) whole_genome.size(); cloc++, rloc++){
				long old = loc_array[cloc];
				if(old == refbase) {
					break;
				} //Already filled with present value
				if(misses > 0 && old >= 0) {
					break;
				} //Already filled with something that has no errors
				char c = toupper(read[cloc]);
				char r = toupper(whole_genome[rloc]);

				if(c == r){
					if(old < 0 || refbase == center_loc){ //If the cell is empty or this key corresponds to center
						loc_array[cloc]=refbase;
					}
				}else{
					misses++;
					if(old >= 0) {
						break;
					} //Already filled with something that has no errors
				}
			}
		}
	}

	//Change 'N' to -2, but only for nocalls, not norefs.  Much faster.
	char nb = 'N';
	for(unsigned int i = 0; i < read.size(); i++){
		if(read[i] == nb) {
			loc_array[i] = -2;
		}
	}
//
	//int minsite = loc_array[0];
	//int maxsite = loc_array[loc_array.size()-1];
	long maxsite = 0;
	for(unsigned int i = loc_array.size()-1; i > 0; i--) {
		if(loc_array[i] > -1) {
			maxsite = loc_array[i];
			break;
		}
	}


	long minsite = 0;
	for(unsigned int i = 0; i < loc_array.size(); i++) {
		if(loc_array[i] > -1) {
			minsite = loc_array[i];
			break;
		}
	}

	//cout << "minsite: " << minsite << endl;
	//cout << "maxsite: " << maxsite << endl;
	for(unsigned int i = 0; i < loc_array.size(); i++) {
		if(loc_array[i] > -1 && (loc_array[i] < minsite || loc_array[i] > maxsite)) {
			loc_array[i] = -1;
		}
	}
/*
	for(unsigned int i = 0; i < loc_array.size(); i++) {
		out << loc_array[i] << " ";
	}
	out << endl;
*/

	int score = calcAffineScore(config, loc_array, read);
	//cout << "score: " << score << endl;
	return score;
}

int maxImperfectScore(Config &config, string &read){
	int maxQ=calcMaxScore(config, read);
	int maxI=maxQ+min(config.POINTS_DEL, config.POINTS_INS-config.POINTS_MATCH2);
	assert(maxI<(maxQ-(config.POINTS_MATCH2+config.POINTS_MATCH2)+(config.POINTS_MATCH+config.POINTS_SUB)));
	return maxI;
}
