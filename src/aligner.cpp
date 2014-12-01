#include "config.h"
#include "genome_index.h"
#include "heap.h"

struct Result {
	int start;
	int stop;
	string matchString;
	Result(int start_, int stop_) {
		start = start_;
		stop = stop_;
	}
};

struct SiteScore {
	int start;
	int stop;
	int score;
	int hits;
	bool perfect;
	int strand;
	vector<int> gapArray;
	SiteScore() {
		start = -1;
		stop = -1;
		score = -1;
		hits = -1;
		perfect = false;
		strand = -1;
	}
	SiteScore(int start_, int stop_, int score_, int hits_, bool perfect_, int strand_, vector<int> gapArray_) {
		start = start_;
		stop = stop_;
		score = score_;
		hits = hits_;
		perfect = perfect_;
		strand = strand_;
		gapArray = gapArray_;
	}
};

bool isNull(SiteScore &ss) {
	if(ss.start == -1 && ss.stop == -1 && ss.score == -1) return true;
	else return false;
}

bool isNull(Triplet &t){
	if(t.column == -1 && t.row == -1 && t.site == -1) return true;
	else return false;
}

void readReads(vector<string> &reads) {
	ifstream ifs(reads_file_name.c_str());
	string line;

	while(getline(ifs, line)) {
		if(line[0] != '>') {
			reads.push_back(line);
		}
	}
}


/*
 *
 * OFFSETS
 *
 */

int getDesiredKeyNumber(int readlen, float density) {
	int slots = readlen - KEYLEN +1;
	int desired = (int) ceil(( readlen * density) / KEYLEN);
	desired = max(minKeysDesired, desired);
	desired = min(slots, desired);
	return desired;
}

void makeOffsets(string read, vector<int> &offsets) {
	float keyDen2 = (( maxDesiredKeys * KEYLEN ) / (float) read.size());
	keyDen2 = max(minKeyDensity, keyDen2);
	keyDen2 = min(keyDensity, keyDen2);
	float keyDen3;
	if(read.size() <= 50){
		keyDen3 = maxKeyDensity;
	}else if(read.size() >= 200){
		keyDen3 = maxKeyDensity-0.5f;
	}else{
		keyDen3 = maxKeyDensity - 0.003333333333f * (read.size()-50);
	}
	keyDen3 = max(keyDensity, keyDen3);

	int desiredKeysNumber = getDesiredKeyNumber(read.size(), keyDen2);

	float interval = (read.size() - KEYLEN) / (float) (max(desiredKeysNumber - 1, 1));
	float f = 0;

	for(int i=0, j = 0; i < desiredKeysNumber; i++){
		if(((unsigned int) round(f+interval)) > read.size()) {
			break;
		}
		offsets.push_back(j);

		f += interval;
		j = min(((int) read.size() - KEYLEN), (max(j+1, (int) round(f))));
	}
}

/*
 *
 * KEY CREATION
 *
 */


int getKeyFromKmer(int start, int stop, string &read) {
	int key = 0;
	for(int i = start; i < stop; i++) {
		int code = getCodeFromBase(read[i]);
		if(code < 0) {
			return -1;
		}
		key = ((key<<2) | code);
	}
	return key;
}

void getReadKeys(string read, vector<int> &offsets, vector<int> &keys) {
	for(unsigned int i = 0; i < offsets.size(); i++) {
		keys.push_back(getKeyFromKmer(offsets[i], offsets[i] + KEYLEN, read));
	}
}


/*
 *
 * KEY COUNTS
 *
 */

int reverseComplementBinary(int kmer, int k){
	int out = 0;
	kmer =~ kmer;
	for(int i = 0; i < k; i++){
		out = ((out<<2) | (kmer&3));
		kmer >>= 2;
	}
	return out;
}

int countKeyHits(int key, int *sizes) {
	int x = sizes[key+1] - sizes[key];
	int rkey = reverseComplementBinary(key, KEYLEN);
	if(key == rkey) return x;
	else {
		int y = sizes[rkey+1] - sizes[rkey];
		return x + y;
	}
}

int countHits(vector<int> &keys, int *sizes, int max_len) {
	vector<int> keys_temp = keys;
	int num_hits = 0;
	for(unsigned int i = 0; i < keys_temp.size(); i++){
		int key = keys_temp[i];
		if(key >= 0){
			int len = countKeyHits(key, sizes);
			if(len > 0 && len < max_len){
				num_hits++;
			}
			else {
				keys[i] = -1;
			}
		}
	}
	return num_hits;
}

int trimHits(vector<int> &keys, vector<int> &keys_temp, int *sizes, int num_hits) {
	if(num_hits > 0) {
		int trigger = (3 * keys_temp.size())/4;
		if(num_hits < 4 && num_hits < trigger) {
			for(unsigned int i = 0; i < keys_temp.size(); i++) {
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, (MAX_LEN * 3) / 2);
		}
		if(num_hits < 3 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, MAX_LEN*2);
		}
		if(num_hits < 3 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, MAX_LEN*3);
		}
		if(num_hits < 2 && num_hits < trigger){
			for(unsigned int i=0; i<keys.size(); i++){
				keys_temp[i] = keys[i];}
			num_hits = countHits(keys_temp, sizes, MAX_LEN*5);
		}
	}
	return num_hits;
}

int getHits(vector<int> &keys, int max_len, vector<int> &starts, vector<int> &stops, int *sizes){
	int num_hits = 0;
	for(unsigned int i = 0; i < keys.size(); i++) {
		int key = keys[i];
		starts.push_back(-1);
		stops.push_back(-1);
		if(key >= 0){
			int len = countKeyHits(key, sizes);
			if(len > 0 && len < max_len){
				starts[i] = sizes[key];
				stops[i] = starts[i] + len;
				num_hits++;
			}
		}
	}
	return num_hits;
}

/*
 *
 * SHRINKS
 *
 */

vector<vector<int> > shrinkArrays(vector<int> &offsets, vector<int> &keys, int num_hits) {

	vector<vector<int> > res;
	vector<int> offsets2;
	vector<int> keys2;

	for(unsigned int i = 0, j = 0; i < keys.size(); i++){
		if(keys[i] >= 0){
			offsets2.push_back(offsets[i]);
			keys2.push_back(keys[i]);
			j++;
		}
	}
	res.push_back(offsets2);
	res.push_back(keys2);
	return res;
}

vector<vector<int> > shrink(vector<int> &starts, vector<int> &stops, vector<int> &offsets, vector<int> &keys, int len){
	int numHits=0;
	for(int i=0; i<len; i++){
		if(starts[i]>=0){numHits++;}
	}
	vector<vector<int> > res;

	if(numHits == (int)offsets.size()){
		return res;
	}else{
		vector<int> offsets2;
		vector<int> starts2;
		vector<int> stops2;

		for(int i=0, j=0; i<len; i++){
			if(starts[i]>=0){
				starts2[j]=starts[i];
				stops2[j]=stops[i];
				offsets2[j]=offsets[i];
				j++;
			}
		}
		res.push_back(starts2);
		res.push_back(stops2);
		res.push_back(offsets2);
		return res;
	}
}

/*
 *
 * ALL BASE COVERED
 *
 */

bool checkIfAllBasesCovered(vector<int> &offsets, string &read) {
	bool all_bases_covered = true;
	if(offsets[0] != 0) {
		all_bases_covered = false;
	}
	else if(offsets[(int)offsets.size() - 1] != ((int)read.size() - KEYLEN)) {
		all_bases_covered = false;
	}
	else{
		for(unsigned int i = 1; i < offsets.size(); i++){
			if(offsets[i] > offsets[i-1] + KEYLEN){
				all_bases_covered = false;
				break;
			}
		}
	}
	return all_bases_covered;
}

bool calculatePretendAllBasesCovered(vector<int> &read_keys,vector<int> &read_keys_final,
				vector<int> &offsets, bool all_bases_covered, string &read) {
	bool pretend_all_bases_covered =
		(all_bases_covered ||
		read_keys_final.size() >= read_keys.size()-4 ||
		(read_keys_final.size() >= 9 &&
				(offsets[offsets.size()-1]-offsets[0]+KEYLEN)>max(40, (int)(read.size()*.75f))));
	return pretend_all_bases_covered;
}


bool isFullyDefined(string &read) {
	for(unsigned int i = 0; i < read.size(); i++) {
		if(getCodeFromBase(read[i]) < 0) return false;
	}
	return true;
}

/*
 *
 * SCORE FUNCTIONS
 *
 */

int calcMaxScore(string &read){
	return POINTS_MATCH + (read.size()-1) * POINTS_MATCH2;
}

int calcMaxScoreZ(vector<int> &offsets){
	int score = 0;
	int a0 = -1;
	int b0 = -1;

	for(unsigned int i = 0; i < offsets.size(); i++){
		int a=offsets[i];
		if(b0 < a){
			score += b0-a0;
			a0 = a;
		}
		b0 = a + KEYLEN;
	}
	score += b0 - a0;
	return score * Z_SCORE_MULT;
}

int scoreZ2(vector<int> &locs, int centerIndex, vector<int> &offsets, int num_approx_hits, int num_hits){

	if(num_approx_hits == 1){
		return SCOREZ_1KEY;
	}

	int center = locs[centerIndex];
	int maxLoc = center + MAX_INDEL2;
	int minLoc = max(0, center - MAX_INDEL);
	int score=0;
	int a0 = -1;
	int b0 = -1;

	for(int i = 0; i < num_hits; i++){
		int loc = locs[i];
		if(loc >= minLoc && loc <= maxLoc){
			int a = offsets[i];
			if(b0 < a){
				score += b0 - a0;
				a0 = a;
			}
			b0 = a + KEYLEN;
		}
	}
	score += b0 - a0;
	score = score * Z_SCORE_MULT;
	return score;
}

int scoreY(vector<int> &locs, int center_index, vector<int > &offsets) {
		int center = locs[center_index];

		int rightIndex = -1;
		for(int i = offsets.size() - 1; rightIndex < center_index; i--){
			if(locs[i] == center){
				rightIndex = i;
			}
		}
		return offsets[rightIndex] - offsets[center_index];
}

int scoreLeft(vector<int> &locs, int center_index) {
	int score=0;
	int prev;
	int loc = locs[center_index];

	for(int i = center_index-1; i >= 0; i--){

		if(locs[i] >= 0){
			prev = loc;
			loc = locs[i];

			int offset = abs(loc - prev);
			if(offset <= MAX_INDEL){
				score += BASE_SCORE;
				if(offset != 0){
					int penalty = min(INDEL_PENALTY+INDEL_PENALTY_MULT*offset, MAX_PENALTY_FOR_MISALIGNED_HIT);
					score-=penalty;
				}
			}else{
				loc=prev;
			}
		}
	}
	return score;
}

int scoreRight(vector<int> &locs, int center_index, int num_hits) {
	int score=0;
	int prev;
	int loc = locs[center_index];

	for(int i = center_index+1; i < num_hits; i++){

		if(locs[i] >= 0){
			prev = loc;
			loc = locs[i];
			int offset = abs(loc - prev);

			if(offset <= MAX_INDEL){
				score += BASE_SCORE;

				if(offset != 0){
					int penalty = min(INDEL_PENALTY + INDEL_PENALTY_MULT * offset, MAX_PENALTY_FOR_MISALIGNED_HIT);
					score -= penalty;
				}
			}else{
				loc = prev;
			}
		}
	}
	return score;
}

int quickScore(vector<int> &locs, int center_index, vector<int> &offsets, int num_approx_hits, int num_hits) {

	if(num_approx_hits == 1) {
		return BASE_SCORE;
	}

	int x = BASE_SCORE + scoreLeft(locs, center_index)+
		scoreRight(locs, center_index, num_hits) - center_index;

	int y = Y_SCORE_MULT * scoreY(locs, center_index, offsets);
	return x+y;
}

int calcMaxQuickScore(vector<int> &offsets, vector<int> keys) {
	int x = keys.size() * BASE_SCORE;
	int y = Y_SCORE_MULT * (offsets[offsets.size()-1] - offsets[0]);
	x += calcMaxScoreZ(offsets);
	return x+y;
}

int calcAffineScore(vector<int> &loc_array, string &read) {
	int score = 0;
	int lastLoc = -3;
	int lastValue = -1;
	int timeInMode = 0;

	for(unsigned int i = 0; i < loc_array.size(); i++){
		int loc = loc_array[i];

		if(loc > 0) {//match
			if(loc == lastValue){//contiguous match
				score += (POINTS_MATCH2);
			}else if(loc == lastLoc || lastLoc < 0) {//match after a sub, or first match
				score += (POINTS_MATCH);
			}else if(loc < lastLoc) {//deletion
				score += (POINTS_MATCH);
				score += POINTS_DEL;
				int dif = lastLoc-loc+1;
				if(dif > MINGAP){
					int rem = dif % GAPLEN;
					int div = (dif - GAPBUFFER2) / GAPLEN;
					score += (div * POINTS_GAP);
					dif = rem + GAPBUFFER2;
				}
				if(dif>LIMIT_FOR_COST_5){
					score+=((dif-LIMIT_FOR_COST_5+MASK5)/TIMESLIP)*POINTS_DEL5;
					dif=LIMIT_FOR_COST_5;
				}
				if(dif>LIMIT_FOR_COST_4){
					score+=(dif-LIMIT_FOR_COST_4)*POINTS_DEL4;
					dif=LIMIT_FOR_COST_4;
				}
				if(dif>LIMIT_FOR_COST_3){
					score+=(dif-LIMIT_FOR_COST_3)*POINTS_DEL3;
					dif=LIMIT_FOR_COST_3;
				}
				if(dif>1){
					score+=(dif-1)*POINTS_DEL2;
				}
				timeInMode=1;
			}else if(loc > lastLoc) {//insertion
				score += (POINTS_MATCH + POINTS_INS_ARRAY_C[min(loc-lastLoc, 5)]);

				timeInMode=1;
			}
			lastLoc=loc;
		}else if(loc==-1){//substitution
			if(lastValue < 0 && timeInMode > 0){//contiguous
				timeInMode++;
				int temp = POINTS_SUB_ARRAY[timeInMode];
				score += temp;
			}else{
				score += POINTS_SUB;
				timeInMode=1;
			}
		}else{
			timeInMode=0;
			score += POINTS_NOCALL;
		}
		lastValue=loc;
	}
	return score;
}

int extendScore(string &read, vector<int> &offsets, vector<int> &values, int center_index, vector<int> &loc_array, int num_hits, int num_approx_hits) {

	int center_loc = values[center_index];
	int min_loc = center_loc - MAX_INDEL;
	int max_loc = center_loc + MAX_INDEL2;

	loc_array.clear();
	for(unsigned int i = 0; i < read.size(); i++) {
		loc_array.push_back(-1);
	}

	//First fill in reverse
	for(int i = 0, keynum = 0; i < num_hits; i++){
		int value=values[i];
		if(value >= min_loc && value <= max_loc){
			int refbase = value;
			keynum++;
			int callbase = offsets[i];
			int misses=0;

			for(int cloc = callbase+KEYLEN-1, rloc = refbase+cloc; cloc >= 0 && rloc >= 0 && rloc < (int) whole_genome.size(); cloc--, rloc--){
				int old = loc_array[cloc];
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
		int value = values[i];
		if(value >= min_loc && value <= max_loc){
			int refbase = value;
			int callbase = offsets[i];
			int misses=0;
			for(int cloc = callbase+KEYLEN, rloc = refbase+cloc; cloc < (int) read.size() && rloc < (int) whole_genome.size(); cloc++, rloc++){
				int old = loc_array[cloc];
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

	int score = calcAffineScore(loc_array, read);
	return score;
}

/*
 *
 * HELP FUNCTIONS
 *
 */

int calcApproxHitsCutoff(int keys, int max_hits, int current_cutoff, bool perfect) {

	int reduction = min( max( (max_hits)/HIT_REDUCTION_DIV, MAX_HITS_REDUCTION2), max(MAXIMUM_MAX_HITS_REDUCTION, keys/8));
	int approxCuttof = max(max(MIN_APPROX_HITS_TO_KEEP, current_cutoff), max_hits-reduction);
	return approxCuttof;
}

bool overlap(int a1, int b1, int a2, int b2){
	return a2<=b1 && b2>=a1;
}

void makeGapArray(vector<int> &locArray, int minLoc, int minGap, vector<int> &gapArray){
	int gaps = 0;
	bool doSort=false;

	if(locArray[0] < 0) { locArray[0] = minLoc; }
	for(unsigned int i = 1; i < locArray.size(); i++){
		if(locArray[i] < 0){ locArray[i] = locArray[i-1] + 1; }
		else { locArray[i] += i; }
		if(locArray[i] < locArray[i-1]) { doSort = true; }
	}

	if(doSort){
		sort(locArray.begin(), locArray.begin() + locArray.size());
	}

	for(unsigned int i = 1; i < locArray.size(); i++){
		int dif = locArray[i] - locArray[i-1];
		if(dif>minGap){
			gaps++;
		}
	}
	if(gaps < 1) {
		locArray.clear();
		return;
	}

	gapArray.push_back(locArray[0]);

	for(unsigned int i = 1, j = 1; i < locArray.size(); i++){
		int dif = locArray[i] - locArray[i-1];
		if(dif > minGap) {
			gapArray.push_back(locArray[i-1]);
			gapArray.push_back(locArray[i]);
			j+=2;
		}
	}

	gapArray.push_back(locArray[locArray.size()-1]);
}


/*
 *
 * PRESCAN FUNCTIONS
 *
 */

void findMaxQscore2(vector<int> &starts, vector<int> &stops, vector<int> &offsets,
		vector<Triplet> &triples, vector<int> &values, vector<int> &keys,
		  int prevMaxHits, bool earlyExit, bool perfectOnly, int *sizes, int *sites, vector<int> &temp){

	int numHits = offsets.size();

	initHeap(offsets.size());
	clear();
	for(int i = 0; i < numHits; i++){
		int start=starts[i];
		int a=sites[start];
		int a2 = a-offsets[i];

		Triplet t(i, start, a2);
		values.push_back(a2);

		add(t);
	}

	int maxQuickScore_=calcMaxQuickScore(offsets, keys);

	int topQscore=-999999999;

	int maxHits=0;

	int approxHitsCutoff;
	int indelCutoff;
	if(perfectOnly){
		approxHitsCutoff=numHits;
		indelCutoff=0;
	}else{
		approxHitsCutoff=max(prevMaxHits, min(MIN_APPROX_HITS_TO_KEEP, numHits-1)); //Faster, same accuracy
		indelCutoff=MAX_INDEL2;
	}

	while(!isEmpty()){
		Triplet t=peek();
		int site=t.site;
		int centerIndex=t.column;

		int maxNearbySite=site;

		int approxHits=0;

		int minsite=site-min(MAX_INDEL, indelCutoff);
		int maxsite=site+MAX_INDEL2;
		for(int column=0, chances=numHits-approxHitsCutoff; column<numHits && chances>=0; column++){
			int x=values[column];

			if(x>=minsite && x<=maxsite){
				if(x > maxNearbySite) {
					maxNearbySite = x;
				}
				approxHits++;
			}else{chances--;}
		}

		if(approxHits>=approxHitsCutoff){

			int qscore=quickScore(values, centerIndex, offsets, approxHits, numHits);

			int scoreZ = scoreZ2(values, centerIndex, offsets, approxHits, numHits);
			qscore+=scoreZ;

			if(qscore>topQscore){

				maxHits=max(approxHits, maxHits);
				approxHitsCutoff=max(approxHitsCutoff, approxHits-1); //Best setting for pre-scan

				topQscore=qscore;

				if(qscore>=maxQuickScore_){
					if(earlyExit){
						temp.push_back(topQscore);
						temp.push_back(maxHits);
						return;
					}
				}
			}
		}

		while(peek().site==site){ //Remove all identical elements, and add subsequent elements
			Triplet t2=poll();

			int row=t2.row+1, col=t2.column;
			if(row<stops[col]){
				t2.row=row;
				int a=sites[row];
				int a2=a-offsets[col];

				t2.site=a2;
				values[col]=a2;
				add(t2);
			}else if(earlyExit && (perfectOnly || size()<approxHitsCutoff)){
				temp.push_back(topQscore);
				temp.push_back(maxHits);
				return;
			}
			if(isEmpty()){break;}
		}

	}
	temp.push_back(topQscore);
	temp.push_back(maxHits);
	return;
}


vector<vector<int> > prescanAllBlocks(int bestScores[], vector<int> &keysP, vector<int> &offsetsP,
		vector<int> &keysM, vector<int> &offsetsM, bool allBasesCovered, int *sizes, int *sites){

//	cout << "offsets.size(): " << offsetsP.size() << ", offsetsM.size(): " << offsetsM.size() << "\n";

	int bestqscore = 0;
	int maxHits = 0;
	int minHitsToScore = MIN_APPROX_HITS_TO_KEEP;

	int maxQuickScore_=calcMaxQuickScore(offsetsP, keysP);

	vector<int> counts;
	vector<int> scores;
	vector<vector<int> > ret;

	vector<int> keys = keysP;
	vector<int> offsets = offsetsP;

	vector<int> starts;
	vector<int> stops;

	int numHits = getHits(keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
		scores.push_back(-9999);
		counts.push_back(0);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<int> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				offsets = r[2];
			}
		}

		clear();
		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;
		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		scores.push_back(temp[0]);
		counts.push_back(temp[1]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			ret.push_back(counts);
			ret.push_back(scores);
			return ret;
		}
	}
	keys.clear();
	keys = keysM;
	offsets.clear();
	offsets = offsetsM;

	starts.clear();
	stops.clear();

	numHits = getHits(keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
		scores.push_back(-9999);
		counts.push_back(0);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<int> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				offsets = r[2];
			}
		}

		clear();
		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;
		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		scores.push_back(temp[0]);
		counts.push_back(temp[1]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			ret.push_back(counts);
			ret.push_back(scores);
			return ret;
		}
	}


	bestScores[1]=max(bestScores[1], maxHits);
	bestScores[3]=max(bestScores[3], bestqscore);

	ret.push_back(counts);
	ret.push_back(scores);
	return ret;
}

/*
 *
 * REVERSE FUNCTIONS
 *
 */

void reverseComplementKeys(vector<int> &reversed_keys, vector<int> &keys ,int k){
	for(unsigned int i = 0, x = keys.size()-1; i < keys.size(); i++){
		int key = reverseComplementBinary(keys[x-i], k);
		reversed_keys.push_back(key);
	}
}

void reverseOffsets(vector<int> &offsets_reversed, vector<int> &offsets, int k, int readlen) {

	for(unsigned int i = 0; i < offsets.size(); i++){
		int x = offsets[offsets.size()-i-1];
		x = readlen - (x+k);
		offsets_reversed.push_back(x);
	}
}

void reverseComplementRead(string &out, string &in){
	for(int i = in.size() - 1; i >= 0; i--){
		if(toupper(in[i]) == 'A') {
			out.push_back('T');
		}
		else if(toupper(in[i]) == 'C') {
			out.push_back('G');
		}
		else if(toupper(in[i]) == 'G') {
			out.push_back('C');
		}
		else {
			out.push_back('A');
		}
	}
}

/*
 *
 * MAKE MATCH STRING
 *
 */

int scoreNoIndelsAndMakeMatchString(string &read, string &ref, int refStart, string &match){
	int score = 0;
	int mode = -1;
	int timeInMode = 0;

	int readStart = 0;
	int readStop = read.size();
	int refStop = refStart + read.size();
	if(refStart < 0 || refStop > (int) ref.size()) {return -99999;}
	if(refStart < 0){
		readStart = 0 - refStart;
		score += POINTS_NOREF * readStart;
	}
	if(refStop > (int) ref.size()){
		int dif = (refStop - ref.size());
		readStop -= dif;
		score += POINTS_NOREF * dif;
	}

	for(int i = readStart; i < readStop; i++){
		char c = read[i];
		char r = ref[refStart+i];

		if(c == r && c != 'N'){
			if(mode == MODE_MS){
				timeInMode++;
				score += POINTS_MATCH2;
			}else{
				timeInMode = 0;
				score += POINTS_MATCH;
			}
			match.push_back('m');
			mode = MODE_MS;
		}else if(c<0 || c=='N'){
			score+=POINTS_NOCALL;
			match.push_back('N');
		}else if(r<0 || r=='N'){
			score+=POINTS_NOREF;
			match.push_back('N');
		}else{
			match.push_back('S');
			if(mode == MODE_SUB) {timeInMode++;}
			else{timeInMode=0;}

			score+=(POINTS_SUB_ARRAY[timeInMode+1]);

			mode=MODE_SUB;
		}
	}
	return score;
}

int makeGref(string &ref, vector<int> &gaps, int refStartLoc, int refEndLoc, string &gref) {

	int g0_old = gaps[0];
	int gN_old = gaps[gaps.size()-1];
	gaps[0] = min(gaps[0], refStartLoc);
	gaps[gaps.size()-1] = max(gN_old, refEndLoc);

	int gpos = 0;

	for(unsigned int i=0; i<gaps.size(); i+=2){
		int x = gaps[i];
		int y = gaps[i+1];

		for(int r=x; r<=y; r++, gpos++){
			gref.push_back(ref[r]);
		}

		if(i+2 < gaps.size()){
			int z = gaps[i+2];
			int gap = z-y-1;

			int rem = gap%GAPLEN;
			int lim = y+ GAPBUFFER +rem;

			int div = (gap-GAPBUFFER2)/GAPLEN;

			for(int r=y+1; r<=lim; r++, gpos++){
				gref.push_back(ref[r]);
			}
			for(int g=0; g<div; g++, gpos++){
				gref.push_back(GAPC);
			}
			for(int r=z-GAPBUFFER; r<z; r++, gpos++){
				gref.push_back(ref[r]);
			}
		}
	}
	int greflimit = gpos;

	gaps[0] = g0_old;
	gaps[gaps.size()-1]=gN_old;

	return greflimit;
}

void fillUnlimited(string &read, string &ref, int refStartLoc, int refEndLoc, int max[], vector<vector<vector<int> > > &packed) {
	int rows = read.size();
	int columns = refEndLoc-refStartLoc+1;

	int maxGain = (read.size()-1) * POINTSoff_MATCH2+POINTSoff_MATCH;
	int subfloor = 0-2*maxGain;
	int BARRIER_I2 = rows-BARRIER_I1;
	int BARRIER_I2b = columns-1;
	int BARRIER_D2 = rows-BARRIER_D1;

	for(int row=1; row<=rows; row++){
		for(int col=1; col<=columns; col++){
			char call0=(row<2 ? '?' : toupper(read[row-2]));
			char call1=toupper(read[row-1]);
			char ref0=(col<2 ? '!' : toupper(ref[refStartLoc+col-2]));
			char ref1=toupper(ref[refStartLoc+col-1]);

			bool match=(call1==ref1 && ref1!='N');
			bool prevMatch=(call0==ref0 && ref0!='N');

			bool gap=(ref1==GAPC);

			if(gap){
				packed[MODE_MS][row][col]=subfloor;
			}else{//Calculate match and sub scores

				int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				int streak=(packed[MODE_MS][row-1][col-1]&TIMEMASK);

				{//Calculate match/sub score

					if(match){

						int scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);
						int scoreD=scoreFromDel+POINTSoff_MATCH;
						int scoreI=scoreFromIns+POINTSoff_MATCH;

						int score;
						int time;
						if(scoreMS>=scoreD && scoreMS>=scoreI){
							score=scoreMS;
							time=(prevMatch ? streak+1 : 1);
						}else if(scoreD>=scoreI){
							score=scoreD;
							time=1;
						}else{
							score=scoreI;
							time=1;
						}

						if(time>MAX_TIME){time=MAX_TIME-MASK5;}
						packed[MODE_MS][row][col]=(score|time);

					}else{

						int scoreMS;
						if(ref1!='N' && call1!='N'){
							scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
								POINTSoff_SUB_ARRAY[streak+1]);
						}else{
							scoreMS=scoreFromDiag+POINTSoff_NOCALL;
						}

						int scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
						int scoreI=scoreFromIns+POINTSoff_SUB;

						int score;
						int time;
						if(scoreMS>=scoreD && scoreMS>=scoreI){
							score=scoreMS;
							time=(prevMatch ? 1 : streak+1);
						}else if(scoreD>=scoreI){
							score=scoreD;
							time=1;
						}else{
							score=scoreI;
							time=1;
						}

						if(time>MAX_TIME){time=MAX_TIME-MASK5;}

						packed[MODE_MS][row][col]=(score|time);

					}
				}
			}

			if(row<BARRIER_D1 || row>BARRIER_D2){
				packed[MODE_DEL][row][col]=subfloor;
			}else{//Calculate DEL score

				int streak=packed[MODE_DEL][row][col-1]&TIMEMASK;

				int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
				int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;

				int scoreMS=scoreFromDiag+POINTSoff_DEL;
				int scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL :
					streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
						streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
							streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
								((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));

				if(ref1=='N'){
					scoreMS+=POINTSoff_DEL_REF_N;
					scoreD+=POINTSoff_DEL_REF_N;
				}else if(gap){
					scoreMS+=POINTSoff_GAP;
					scoreD+=POINTSoff_GAP;
				}

				int score;
				int time;
				if(scoreMS>=scoreD){
					score=scoreMS;
					time=1;
				}else{
					score=scoreD;
					time=streak+1;
				}

				if(time>MAX_TIME){time=MAX_TIME-MASK5;}
				packed[MODE_DEL][row][col]=(score|time);
			}

			//Calculate INS score
			if(gap || (row<BARRIER_I1 && col>1) || (row>BARRIER_I2 && col<BARRIER_I2b)){
				packed[MODE_INS][row][col]=subfloor;
			}else{//Calculate INS score

				int streak=packed[MODE_INS][row-1][col]&TIMEMASK;

				int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
				int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;

				int scoreMS=scoreFromDiag+POINTSoff_INS;
				int scoreI=scoreFromIns+POINTSoff_INS_ARRAY[streak+1];

				int score;
				int time;
				if(scoreMS>=scoreI){
					score=scoreMS;
					time=1;
				}else{
					score=scoreI;
					time=streak+1;
				}

				if(time>MAX_TIME){time=MAX_TIME-MASK5;}
				packed[MODE_INS][row][col]=(score|time);
			}
		}
	}

	int maxCol=-1;
	int maxState=-1;
	int maxScore = std::numeric_limits<int>::min();

	for(int state=0; state<3; state++){
		for(int col=1; col<=columns; col++){
			int x=packed[state][rows][col]&SCOREMASK;
			if(x>maxScore){
				maxScore=x;
				maxCol=col;
				maxState=state;
			}
		}
	}
	maxScore>>=SCOREOFFSET;

	max[0] = rows;
	max[1] = maxCol;
	max[2] = maxState;
	max[3] = maxScore;

}

void traceback(string &read, string &ref, int refStartLoc, int refEndLoc, int row, int col, int state, string &out, vector<vector<vector<int> > > &packed) {

	int outPos=0;
	int gaps=0;

	//string out;

	while(row>0 && col>0){

		//cout << out.capacity() << " ";
		//cout << out.size() << " ";
		//cout << col << " " << endl;
		/*if(col == 143) {
			//cout << out << endl;
			cout << out.capacity() << endl;
			cout << out.size() << endl;
			if(out.capacity() == out.size()) {
				out.reserve(out.size()+1);
			}
		}*/

		int time=packed[state][row][col]&TIMEMASK;
		int prev;

		if(state==MODE_MS){
			if(time>1){prev=(int)state;}
			else{
				int scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				int scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				int scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
				else{prev=MODE_INS;}
			}

			char c=toupper(read[row-1]);
			char r=toupper(ref[refStartLoc+col-1]);
			if(c==r){
				out.push_back('m');
				//cout << 'm';
			}else{
				out.push_back('s');
				//cout << 'S';
			}

			row--;
			col--;
		}else if(state==MODE_DEL){
			if(time>1){prev=(int)state;}
			else{
				int scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
				int scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
				else{prev=MODE_DEL;}
			}

			char r=toupper(ref[refStartLoc+col-1]);
			if(r==GAPC){
				out.push_back('-');
				//cout << '-';
				gaps++;
			}else{
				out.push_back('d');
				//cout << 'D';
			}
			col--;
		}else{
			if(time>1){prev=(int)state;}
			else{
				int scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
				int scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
				if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else{prev=MODE_INS;}
			}

			if(col==0){
				out.push_back('x');
				//cout << 'X';
			//}else if(col>=5000){
		//		out.push_back('Y'); //TODO
			}else{
				out.push_back('i');
				//cout << 'I';
			}
			row--;
		}

		state=prev;
		outPos++;
	}

	if(col!=row){
		while(row>0){
			out.push_back('x');
			outPos++;
			row--;
			col--;
		}
		if(col>0){
			//do nothing
		}
	}

	reverse(out.begin(),out.end());

	if(gaps==0){
		/*for(unsigned int j = 0; j < out.size(); j++) {
			cout << out[j];
		}
		cout << "\n";
*/
		return;
	}

	string out3;
	int j = 0;
	for(unsigned int i=0; i<out.size(); i++){
		char c=out[i];
		if(c!=GAPC){
			out3.push_back(c);
			j++;
		}else{
			int lim=j+GAPLEN;
			for(; j<lim; j++){
				out3.push_back('d');
			}
		}
	}
	out = out3;
/*
	for(unsigned int j = 0; j < out3.size(); j++) {
		cout << out3[j];
	}
	cout << "\n";
*/
	return;
}


void fillLimited(string &read, string &ref, int refStartLoc, int refEndLoc, int score, vector<int> &gaps, int max[], Result &r) {
	string gref;
	int grefLimit = (refEndLoc - refStartLoc) + 2 * MARIC_PADDING;
	if((int)gaps.size() == 0) {
		for(int i = refStartLoc-MARIC_PADDING; i < refEndLoc + MARIC_PADDING; i++) {
			gref.push_back(whole_genome[i]);
		}
	}
	else {
		grefLimit = makeGref(ref, gaps, (refStartLoc - MARIC_PADDING), (refEndLoc + MARIC_PADDING), gref);
	}
	int rows = read.size();
	int columns = grefLimit+1;
	vector<vector<vector<int> > > packed;
	/*int ***packed;
	packed = new int**[3];
	for(int j = 0; j < 3; j++) {
		packed[j] = new int*[rows+1];
		for(int i = 0; i < rows+1; ++i) {
			packed[j][i] = new int[columns+1];
		}
	}*/

	for(int matrix=0; matrix<3; matrix++){
		vector<int> row;
		for(int i = 0; i <= columns; i++) {
			row.push_back(0);
		}
		vector<vector<int> > mat;
		mat.push_back(row);
		for(int i=1; i<=rows; i++){
			vector<int> row;
			for(int j=0; j<columns+1; j++){
				//packed[matrix][i][j]=BADoff;
				row.push_back(BADoff);
			}
			mat.push_back(row);
		}
		for(int i=0; i<=rows; i++){

			int prevScore=(i<2 ? 0 : mat[i-1][0]);
			int score=prevScore+POINTSoff_INS_ARRAY[i];

			mat[i][0]=score;
		}
		packed.push_back(mat);
	}
/*
	for(int matrix=0; matrix<3; matrix++){
		vector<int> row;
		for(int i = 0; i <= columns; i++) {
			packed[matrix][0][i] = 0;
		}
		for(int i=1; i<=rows; i++){
			for(int j=0; j<columns+1; j++){
				packed[matrix][i][j]=BADoff;

			}

		}
		for(int i=0; i<=rows; i++){

			int prevScore=(i<2 ? 0 : packed[matrix][i-1][0]);
			int score=prevScore+POINTSoff_INS_ARRAY[i];

			packed[matrix][i][0]=score;
		}
	}
*/
	fillUnlimited(read, gref, 0, grefLimit, max, packed);

	//string out;
	//cout << gref << endl;
	/*for(int i = 0; i < rows+1; i++) {
		for(int j = 0; j < columns+1; j++) {
			cout << packed[0][i][j] << " ";
		}
		cout << endl;
	}*/
	//cout << endl;
	traceback(read, gref, 0, grefLimit, max[0], max[1], max[2], r.matchString, packed);

/*
	for(int i = 0; i < rows+1; i++) {
		for(int j = 0; j < columns+1; j++) {
			cout << packed[0][i][j] << " ";
		}
		cout << endl;
	}
*/
/*
	for(int j = 0; j < 3; j++) {
		for(int i = 0; i < rows+1; i++) {
			delete [] packed[j][i];
		}
		delete [] packed[j];
	}
	delete [] packed;
	*/
}

void makeMatchStringForSite(SiteScore ss, string &read, int *sizes, int *sites, Result &r) {
	string match;
	int score = scoreNoIndelsAndMakeMatchString(read, whole_genome, ss.start, match);

	int max[4];
	//fillLimited(read, whole_genome, ss.start, ss.stop, score, ss.gapArray, max, r);
}

void makeMatchString(vector<SiteScore> &results, string &read, string &read_reverse, int *sizes, int *sites, vector<Result> &resultsFinal) {
	int max_score = -9999;
	SiteScore ss;
	for(unsigned int i = 0; i < results.size(); i++) {
		if(results[i].score > max_score) {
			max_score = results[i].score;
			ss = results[i];
		}
	}
	Result r = Result(ss.start, ss.stop);
	if(ss.strand == 0) {
		makeMatchStringForSite(ss, read, sizes, sites, r);
	}
	else {
		makeMatchStringForSite(ss, read_reverse, sizes, sites, r);
	}
	resultsFinal.push_back(r);
}


/*
 *
 * ALIGNER
 *
 */

void align(int bestScores[], vector<int> &keys, string &read, vector<int> &offsets, int *sizes,
		int *sites, vector<SiteScore> &results, bool all_bases_covered,
		int max_score, int max_quick_score, bool fully_defined, int strand) {
	vector<int> starts;
	vector<int> stops;
	int get_hits = getHits(keys, MAX_LEN, starts, stops, sizes);

	bool filter_by_qscore = (get_hits >= 5);
	int minScore = (int) (MIN_SCORE_MULT * max_score);
	int minQuickScore = (int)(MIN_QSCORE_MULT * max_quick_score);

	int current_top_score = bestScores[0];
	int bestqscore = bestScores[3];
	int max_hits = bestScores[1];
	int perfects_found = bestScores[5];

	int cutoff = max(minScore, (int)(current_top_score * DYNAMIC_SCORE_THRESH));
	int qcutoff = max(bestScores[2], minQuickScore);

	int approx_hits_cutoff = calcApproxHitsCutoff(get_hits, max_hits, MIN_APPROX_HITS_TO_KEEP, current_top_score >= max_score);
	if(approx_hits_cutoff > get_hits){
		return;
	}
	bool short_circuit = (all_bases_covered && get_hits == (int)offsets.size() && filter_by_qscore);

	if(current_top_score >= max_score){
		qcutoff = max(qcutoff, (int)(max_quick_score*DYNAMIC_QSCORE_THRESH_PERFECT));
	}

	vector<int> sites_tmp;
	clear();
	initHeap(offsets.size());

	for(int i = 0; i < get_hits; i++){
		int start = starts[i];

		int a = sites[start];
		int a2 = a - offsets[i];

		Triplet t(i, start, a2);
		sites_tmp.push_back(a2);

		add(t);
	}
	if(write_output) cout << "\n";

	vector<int> loc_array;
	SiteScore prevSS;

	int numberOfSites = 0;

	while(!isEmpty()) {
		Triplet t = peek();
		int site = t.site;
		int center_index = t.column;
		int max_nearby_site = site;
		int approx_hits = 0;
		int minsite = site - MAX_INDEL;
		int maxsite = site + MAX_INDEL2;
		for(int column = 0, chances = get_hits - approx_hits_cutoff; column < get_hits && chances >= 0; column++) {
			int x = sites_tmp[column];
			if(x >= minsite && x <= maxsite) {
				if(x > max_nearby_site) {
					max_nearby_site = x;
				}
				approx_hits++;
			} else {
				chances--;
			}
		}

		if(approx_hits >= approx_hits_cutoff){

			int score;
			int qscore = quickScore(sites_tmp, center_index, offsets, approx_hits, get_hits);
			int scoreZ = scoreZ2(sites_tmp, center_index, offsets, approx_hits, get_hits);
			qscore+=scoreZ;

			int mapStart = site;
			int mapStop = max_nearby_site;

			if(qscore < qcutoff){
				score=-1;
			}else{
				if(short_circuit && qscore == max_quick_score){
					score = max_score;
				}else{
					score = extendScore(read, offsets, sites_tmp, center_index, loc_array, get_hits, approx_hits);

					int min = std::numeric_limits<int>::max();
					int max = std::numeric_limits<int>::min();
					for(unsigned int i = 0; i < loc_array.size(); i++){
						int x = loc_array[i];
						if(x >- 1){
							if(x < min) {
								min = x;
							}
							if(x > max) {
								max = x;
							}
						}
					}
					if(min<0 || max<0){
						score=-99999;
					}
					mapStart = min;
					mapStop = max;
				}
				if(score == max_score){
					qcutoff = max(qcutoff, (int)(max_quick_score * DYNAMIC_QSCORE_THRESH_PERFECT));
					approx_hits_cutoff = calcApproxHitsCutoff(get_hits, max_hits, MIN_APPROX_HITS_TO_KEEP, current_top_score >= max_score);
				}

				if(score >= cutoff){
					qcutoff = max(qcutoff, (int)(qscore * DYNAMIC_QSCORE_THRESH));
					bestqscore = max(qscore, bestqscore);
				}
			}
			if(score>=cutoff){

				if(score > current_top_score){
					max_hits = max(approx_hits, max_hits);
					approx_hits_cutoff = calcApproxHitsCutoff(get_hits, max_hits, approx_hits_cutoff, current_top_score >= max_score);
					cutoff = max(cutoff, (int)(score * DYNAMIC_SCORE_THRESH));

					if(score >= max_score){
						cutoff = max(cutoff, (int)(score*0.95f));
					}
					current_top_score = score;
				}

				int site2 = mapStart;
				int site3 = mapStop + read.size() - 1;

				//cout << "start: " << site2 << ", stop: " << site3 << ", size: " << (site3-site2) << endl;
				//cout << "score: " << score << endl;
				//cout << "strand: " << strand << endl;
				vector<int> gapArray;
				if(site3-site2 >= MINGAP + (int)read.size()){

					makeGapArray(loc_array, site2, MINGAP, gapArray);
					if(gapArray.size() != 0){
						gapArray[0] = min(gapArray[0], site2);
						gapArray[gapArray.size()-1] = max(gapArray[gapArray.size()-1], site3);
					}
				}

				SiteScore ss;
				bool perfect1 = score == max_score && fully_defined;
				bool inbounds = (site2 >= 0 && site3 < (int)whole_genome.size());
				//cout << "perfect1: " << perfect1 << ", max_score: " << max_score << " fully_defined: " << fully_defined << endl;

				if(inbounds && (int) gapArray.size() == 0 && isNull(prevSS) == false
						&& strand == prevSS.strand && overlap(prevSS.start, prevSS.stop, site2, site3)){

					int betterScore = max(score, prevSS.score);
					int minStart = min(prevSS.start, site2);
					int maxStop = max(prevSS.stop, site3);
					bool perfect2 = prevSS.score == max_score && fully_defined;

					bool shortEnough = (maxStop-minStart < 2 * (int)read.size());

					if(prevSS.start == site2 && prevSS.stop == site3){
						results[numberOfSites-1].score = betterScore;
						results[numberOfSites-1].perfect = (prevSS.perfect || perfect1 || perfect2);
						//if(prevSS.perfect){prevSS.semiperfect=true;}

					} else if(shortEnough && prevSS.start == site2){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[numberOfSites-1].stop = site3;
							if(!prevSS.perfect) { perfects_found++; }
							results[numberOfSites-1].perfect = true;
							results[numberOfSites-1].score = betterScore;
						}else{
							results[numberOfSites-1].stop = maxStop;
							// prevSS.setPerfect(bases); TODO add perfect set
						}
						results[numberOfSites-1].score = betterScore;
					}else if(shortEnough && prevSS.stop == site3){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[numberOfSites-1].start = site2;
							if(!prevSS.perfect) { perfects_found++; }
							results[numberOfSites-1].perfect = true;
							results[numberOfSites-1].score = betterScore;
						}else{
							results[numberOfSites-1].start = minStart;
							//prevSS.setPerfect(bases); TODO add perfect set
						}
						results[numberOfSites-1].score = betterScore;
					}else if(shortEnough && (maxStop - minStart <= (int) read.size() + MAX_SUBSUMPTION_LENGTH) && !perfect1 && !perfect2){
						results[numberOfSites-1].start = minStart;
						results[numberOfSites-1].stop = maxStop;
						results[numberOfSites-1].score = betterScore;
						//prevSS.setPerfect(bases); TODO add perfect set
					}else{
						ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
						//cout << "first: start: " << site2+8000 << ", stop: " << site3+8000 << "\n";
						//if(!perfect1){ss.setPerfect(bases);} TODO add perfect set
					}

				}else if(inbounds){
					ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
					//cout << "second: start: " << site2+8000 << ", stop: " << site3+8000 << "\n";
					//if(!perfect1){ss.setPerfect(bases);} TODO add perfect setter
				}

				if(isNull(ss) == false){
					results.push_back(ss);
					numberOfSites++;
					if(ss.perfect){
						if(isNull(prevSS) == true || !prevSS.perfect || overlap(ss.start, ss.stop, prevSS.start, prevSS.stop)){
							perfects_found++;
							if(perfects_found >= 2) {
								clear();
								break;
							}
						}
					}
					prevSS = ss;
				}
			}
		}

		while(peek().site == site){

			Triplet t2 = poll();
			int row = t2.row+1;
			int col = t2.column;
			if(row < stops[col]){
				t2.row=row;
				int a = sites[row];
				int a2 = a-offsets[col];
				t2.site = a2;
				sites_tmp[col]=a2;
				add(t2);
			}
			else if(size() < approx_hits_cutoff){
				return;
			}
		}
	}
}

void processRead(int *sizes, int *sites, string &read, vector<Result> &resultsFinal) {

	string read_reverse;
	reverseComplementRead(read_reverse, read);

	vector<int> offsets;
	makeOffsets(read, offsets);
	vector<int> read_keys;
	getReadKeys(read, offsets, read_keys);
	vector<int> read_keys_copy = read_keys;

	int num_hits = countHits(read_keys_copy, sizes, MAX_LEN);
	num_hits = trimHits(read_keys, read_keys_copy, sizes, num_hits);

	vector<vector<int> > res;

	//TODO add this if if(num_hits < read_keys.size()) {
	res = shrinkArrays(offsets, read_keys_copy, num_hits);

	offsets = res[0];
	vector<int> read_keys_final = res[1];

	vector<SiteScore> results;

	bool all_bases_covered = checkIfAllBasesCovered(offsets, read);
	bool pretend_all_bases_covered = calculatePretendAllBasesCovered(read_keys, read_keys_final, offsets, all_bases_covered, read);
	bool prescan_qscore = (num_hits >= 5);

	int quick_max_score = calcMaxQuickScore(offsets, read_keys_final);

	int bestScores[6] = {0};
	vector<int> precounts;
	vector<int> prescores;

	int hitsCutoff=0;
	int qscoreCutoff = (int)(MIN_QSCORE_MULT*quick_max_score);

	vector<int> offsets_reversed;
	reverseOffsets(offsets_reversed, offsets, KEYLEN, read.size());

	vector<int> keys_reversed;
	reverseComplementKeys(keys_reversed, read_keys_final, KEYLEN);

	if(prescan_qscore){
		vector<vector<int> > prescanResults = prescanAllBlocks(bestScores, read_keys_final, offsets,
				keys_reversed, offsets_reversed, pretend_all_bases_covered, sizes, sites);

		precounts=prescanResults[0];
		prescores=prescanResults[1];

		if(bestScores[1]<MIN_APPROX_HITS_TO_KEEP){return;}
		if(bestScores[3]<quick_max_score*MIN_QSCORE_MULT2){return;} //if(bestScores[3]<maxQuickScore(offsetsP, keyScoresP)*.10f){return result;}

		if(bestScores[3]>=quick_max_score && pretend_all_bases_covered){

			hitsCutoff=calcApproxHitsCutoff(read_keys_final.size(), bestScores[1], MIN_APPROX_HITS_TO_KEEP, true);
			qscoreCutoff=max(qscoreCutoff, (int)(bestScores[3]*DYNAMIC_QSCORE_THRESH_PERFECT));

		}else{
			hitsCutoff=calcApproxHitsCutoff(read_keys_final.size(), bestScores[1], MIN_APPROX_HITS_TO_KEEP, false);
			qscoreCutoff=max(qscoreCutoff, (int)(bestScores[3]*PRESCAN_QSCORE_THRESH));
		}
	}

	int max_score = calcMaxScore(read);

	bool fully_defined = isFullyDefined(read);

	if(precounts[0] >= hitsCutoff && prescores[0] >= qscoreCutoff) {
		align(bestScores, read_keys_final, read, offsets, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 0);
	}
	if(precounts[1] >= hitsCutoff && prescores[1] >= qscoreCutoff) {
		align(bestScores, keys_reversed, read_reverse, offsets_reversed, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 1);
	}

	makeMatchString(results, read, read_reverse, sizes, sites, resultsFinal);

	return;
}

void writeResults(vector<Result> results) {

	ofstream out_res(results_file.c_str());
	int br = 1;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		out_res << br << "-" << r.start << "-" << r.stop << "-" << (r.stop-r.start) << endl;
		out_res << r.matchString << endl;
		br++;
	}
}

int main(int argc, char *argv[]) {
	int **res = readIndex(whole_genome);

	int *sizes = res[1];
	int *sites = res[3];

	vector<string> reads;
	readReads(reads);

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	cout << "NUMBER OF READS: "<< reads.size() << "\n";

	vector<Result> results;
	for(unsigned int i = 0; i < reads.size(); i++) {
		//cout << (i+1) << "------------------- \n";
		processRead(sizes, sites, reads[i], results);
	}
	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "ALIGNED READS: " << timefinal << endl;

	cout << "END" << endl;

	writeResults(results);

	return 0;
}
