#include "config.h"
#include "genome_index.h"
#include "heap.h"

struct Statistic {
	int covered;
	int reads;
	int reads_with_result;
	int gapped;
	int ungapped;
	int gapped_with_start;
	int gapped_with_stop;
	int gapped_with_start_and_stop;
	int gapped_with_precise_start;
	int gapped_with_precise_stop;
	int gapped_with_precise_start_and_stop;
	int ungapped_with_start;
	int ungapped_with_stop;
	int ungapped_with_start_and_stop;
	int ungapped_with_precise_start;
	int ungapped_with_precise_stop;
	int ungapped_with_precise_start_and_stop;
	int exons;
	int correct_exons;
	int same_numbered_exons;
	int correct_exon_start;
	int correct_exon_stop;
	int correct_exon_start_and_stop;
	int percentage;
	int precise_percentage;
	int gapped_percentage;
	int gapped_precise_percentage;
	int ungapped_percentage;
	int ungapped_precise_percentage;
	Statistic(int num_reads) {
		covered = 0,
		reads = num_reads;
		reads_with_result = 0;
		gapped = 0;
		ungapped = 0;
		gapped_with_start = 0;
		gapped_with_stop = 0;
		gapped_with_start_and_stop = 0;
		gapped_with_precise_start = 0;
		gapped_with_precise_stop = 0;
		gapped_with_precise_start_and_stop = 0;
		ungapped_with_start = 0;
		ungapped_with_stop = 0;
		ungapped_with_start_and_stop = 0;
		ungapped_with_precise_start = 0;
		ungapped_with_precise_stop = 0;
		ungapped_with_precise_start_and_stop = 0;
		exons = 0;
		correct_exons = 0;
		same_numbered_exons = 0;
		correct_exon_start = 0;
		correct_exon_stop = 0;
		correct_exon_start_and_stop = 0;
		percentage = 0;
		precise_percentage = 0;
		gapped_percentage = 0;
		gapped_precise_percentage = 0;
		ungapped_percentage = 0;
		ungapped_precise_percentage = 0;
	}
};

struct Result {
	int br;
	int start;
	int stop;
	string matchString;
	vector<int> gapArray;
	int precise_start;
	int precise_stop;
	Result(int br_, int start_, int stop_) {
		br = br_;
		start = start_;
		stop = stop_;
	}
};


struct ThreadData3 {
	int thread_id;
	int *sizes;
	int *sites;
	vector<string> *reads;
	int start;
	int stop;
	string *whole_genome;
	vector<Result> *results;
	ThreadData3() {
	}
	ThreadData3(int thread_id_, int *sizes_, int *sites_, vector<string> *reads_, int start_, int stop_, string *whole_genome_, vector<Result> *results_) {
		thread_id = thread_id_;
		sizes = sizes_;
		sites = sites_;
		reads = reads_;
		start = start_;
		stop = stop_;
		whole_genome = whole_genome_;
		results = results_;
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

void gapArrayFromString(vector<int> &gapArray, string str) {
	char delim = ' ';
	char delim2 = '-';
	stringstream ss(str);
	string item;
	while (std::getline(ss, item, delim)) {
		if(item[item.size()-1] == ',') {
			item.erase(item.size()-1);
		}
		stringstream ss2(item);
		string first;
		string second;
		getline(ss2, first, delim2);
		getline(ss2, second, delim2);
		gapArray.push_back(atoi(first.c_str())-1);
		gapArray.push_back(atoi(second.c_str())-1);
	}
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void readGapArray(Result &r, string line) {
	std::vector<std::string> elems;
	split(line, ' ', elems);
	for(unsigned int i = 0; i < elems.size(); i++) {
		//cout << elems[i] << endl;
		string tmp;
		if(i != elems.size() - 1) {
			tmp = elems[i].substr(0, elems[i].size()-1);
		}
		else {
			tmp = elems[i];
		}
		vector<string> el;
		//cout << tmp << endl;
		split(tmp, '-', el);
		r.gapArray.push_back(atoi(el[0].c_str()));
		r.gapArray.push_back(atoi(el[1].c_str()));
	}
	/*for(unsigned int i = 0; i < r.gapArray.size(); i++) {
		cout << r.gapArray[i] << endl;
	}*/
}

void readReads(vector<string> &reads, map<int, Result> &results, string infile) {
	ifstream ifs(infile.c_str());
	if(!ifs) {
		cout << "File " << infile << " does not exist." << endl;
		exit(-1);
	}
	string line;
	string line2;
	string line3;
	string line4;
	string line5;
	string line6;

	int br = 1;
	while(getline(ifs, line)) {
		getline(ifs, line2);
		getline(ifs, line3);
		getline(ifs, line4);
		getline(ifs, line5);
		getline(ifs, line6);
		//int br = atoi(line.c_str());
		Result r = Result(br, atoi(line2.c_str()), atoi(line3.c_str()));
		//readGapArray(r, line5);
		r.matchString = line4;
		gapArrayFromString(r.gapArray, line5);
		//results.push_back(r);
		std::map<int,Result>::iterator it = results.begin();
		results.insert(it, pair<int,Result>(br, r));
		reads.push_back(line6);
		total_base_num += line6.size();
		br++;
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
	int x;
	if(key >= length_of_sizes) {
		//cout << "length: "  << length_of_sizes << endl;
		//cout << "key: " << key << endl;
	}
	if(key+1 == length_of_sizes) {
		x = length_of_sites - sizes[key];
	}
	else {
		x = sizes[key+1] - sizes[key];
	}
	int rkey = reverseComplementBinary(key, KEYLEN);
	if(key == rkey) return x;
	else {
		int y;
		if(rkey+1 == length_of_sizes) {
			y = length_of_sites - sizes[rkey];
		}
		else {
			y = sizes[rkey+1] - sizes[rkey];
		}
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
		if(num_hits < 20 && num_hits < trigger) {
			for(unsigned int i = 0; i < keys_temp.size(); i++) {
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, (MAX_LEN * 3) / 2);
		}
		if(num_hits < 18 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, MAX_LEN*2);
		}
		if(num_hits < 16 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(keys_temp, sizes, MAX_LEN*3);
		}
		if(num_hits < 14 && num_hits < trigger){
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

int extendScore(string &read, vector<int> &offsets, vector<int> &values, int center_index, vector<int> &loc_array, int num_hits, int num_approx_hits, string &whole_genome) {

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

int maxImperfectScore(string &read){
	int maxQ=calcMaxScore(read);
	int maxI=maxQ+min(POINTS_DEL, POINTS_INS-POINTS_MATCH2);
	assert(maxI<(maxQ-(POINTS_MATCH2+POINTS_MATCH2)+(POINTS_MATCH+POINTS_SUB)));
	return maxI;
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
	Heap heap = Heap();
	heap.initHeap(offsets.size());
	heap.clear();
	for(int i = 0; i < numHits; i++){
		int start=starts[i];
		int a=sites[start];
		int a2 = a-offsets[i];

		Triplet t(i, start, a2);
		values.push_back(a2);

		heap.add(t);
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

	while(!heap.isEmpty()){
		Triplet t=heap.peek();
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

		while(heap.peek().site==site){ //Remove all identical elements, and add subsequent elements
			Triplet t2=heap.poll();

			int row=t2.row+1, col=t2.column;
			if(row<stops[col]){
				t2.row=row;
				int a=sites[row];
				int a2=a-offsets[col];

				t2.site=a2;
				values[col]=a2;
				heap.add(t2);
			}else if(earlyExit && (perfectOnly || heap.size()<approxHitsCutoff)){
				temp.push_back(topQscore);
				temp.push_back(maxHits);
				return;
			}
			if(heap.isEmpty()){break;}
		}

	}
	temp.push_back(topQscore);
	temp.push_back(maxHits);
	return;
}

void prescanAllBlocks(vector<vector<int> > &prescanResults, int bestScores[], vector<int> &keysP, vector<int> &offsetsP,
		vector<int> &keysM, vector<int> &offsetsM, bool allBasesCovered, int *sizes, int *sites){

//	cout << "offsets.size(): " << offsetsP.size() << ", offsetsM.size(): " << offsetsM.size() << "\n";

	int bestqscore = 0;
	int maxHits = 0;
	int minHitsToScore = MIN_APPROX_HITS_TO_KEEP;

	int maxQuickScore_=calcMaxQuickScore(offsetsP, keysP);

	vector<int> counts;
	vector<int> scores;
	prescanResults.push_back(counts);
	prescanResults.push_back(scores);
	//vector<vector<int> > ret;

	vector<int> keys = keysP;
	vector<int> offsets = offsetsP;

	vector<int> starts;
	vector<int> stops;

	int numHits = getHits(keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
		//scores.push_back(-9999);
		//counts.push_back(0);
		prescanResults[0].push_back(0);
		prescanResults[1].push_back(-9999);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<int> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				offsets = r[2];
			}
		}

		//heap.clear();
		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;
		//cout << "findMaxQscore2" << endl;

		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		//scores.push_back(temp[0]);
		//counts.push_back(temp[1]);
		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			//prescanResults.push_back(counts);
			//prescanResults.push_back(scores);
			return;
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
		//scores.push_back(-9999);
		//counts.push_back(0);
		prescanResults[0].push_back(0);
		prescanResults[1].push_back(-9999);
	}else{

		if(numHits < (int) keys.size()){
			vector<vector<int> > r = shrink(starts, stops, offsets, keys, offsets.size());
			if(r.size() != 0){
				starts = r[0];
				stops = r[1];
				offsets = r[2];
			}
		}

		//heap.clear();
		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;
		//cout << "findMaxQscore22" << endl;
		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		//scores.push_back(temp[0]);
		//counts.push_back(temp[1]);

		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			//prescanResults.push_back(counts);
			//prescanResults.push_back(scores);
			return;
		}
	}


	bestScores[1]=max(bestScores[1], maxHits);
	bestScores[3]=max(bestScores[3], bestqscore);

	//prescanResults.push_back(counts);
	//prescanResults.push_back(scores);
	return;
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
			match.push_back('s');
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

void traceback(string &read, string &ref, int refStartLoc, int refEndLoc, int row, int col, int state, Result &r, vector<vector<vector<int> > > &packed) {

	/*if(r.br == 26) {
		cout << "here" << endl;
	}*/

	int outPos=0;
	int gaps=0;

	int addon = refEndLoc - col;

	string out;

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

	//cout << read << endl;
	//cout << ref << endl;

	/*if(r.br == 26) {
		cout << "er" << endl;
	}
*/
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


	r.precise_start = r.start-MARIC_PADDING + col;
	r.precise_stop = r.stop+MARIC_PADDING -addon;

	if(gaps==0){
		/*for(unsigned int j = 0; j < out.size(); j++) {
			cout << out[j];
		}
		cout << "\n";
*/
		r.matchString = out;
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
	r.matchString = out3;
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
	bool gapped = true;
	if((int)gaps.size() == 0) {
		gapped = false;
		for(int i = refStartLoc-MARIC_PADDING; i < refEndLoc + MARIC_PADDING; i++) {
			gref.push_back(ref[i]);
		}
	}
	else {
		/*for(int i = 0; i < gaps.size(); i++) {
			cout << gaps[i] << " ";
		}
		cout << endl;*/
		grefLimit = makeGref(ref, gaps, (refStartLoc - MARIC_PADDING), (refEndLoc + MARIC_PADDING), gref);
	}
	//cout << gref << endl;
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

	ofstream test("tmp3/test.out");
	for(unsigned int i = 0; i < packed[0].size(); i++) {
		for(unsigned int j = 0; j < packed[0][i].size(); j++) {
			test << packed[0][i][j] << " ";
		}
		test << endl;
	}

	//string out;
	//cout << gref << endl;
	/*for(int i = 0; i < rows+1; i++) {
		for(int j = 0; j < columns+1; j++) {
			cout << packed[0][i][j] << " ";
		}
		cout << endl;
	}*/
	//cout << endl;
	traceback(read, gref, 0, grefLimit, max[0], max[1], max[2], r, packed);

	if(r.gapArray.size() == 0) r.precise_stop--;
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

void makeMatchStringForSite(SiteScore ss, string &read, int *sizes, int *sites, Result &r, string &whole_genome, int threadId) {

	if(ss.perfect) {
		//r.matchString = "mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm";
		r.start = ss.start;
		r.stop = ss.stop;
		r.precise_start = r.start;
		r.precise_stop = r.stop;
		for(unsigned int i = 0; i < read.size(); i++) {
			r.matchString.push_back('m');
		}
		return;
	}

	string match;
	int score = scoreNoIndelsAndMakeMatchString(read, whole_genome, ss.start, match);
	int imperfectScore = maxImperfectScore(read);
	if(score > imperfectScore) {
		r.start = ss.start;
		r.stop = ss.stop;
		r.precise_start = r.start;
		r.precise_stop = r.stop;
		r.matchString = match;
		return;
	}

	int max[4];
	fillLimited(read, whole_genome, ss.start, ss.stop, score, ss.gapArray, max, r);
}

void makeMatchString(vector<SiteScore> &results, string &read, string &read_reverse, int *sizes, int *sites, vector<Result> &resultsFinal, string &whole_genome, int threadId, int br) {
	int max_score = -9999;
	SiteScore ss;

	for(unsigned int i = 0; i < results.size(); i++) {
		if(results[i].score > max_score) {
			max_score = results[i].score;
			SiteScore ress = results[i];
			ss = results[i];
		}
		/*if(results[i].stop == 144497458) {
			ss = results[i];
		}*/
	}
	Result r = Result(br, ss.start, ss.stop);
	r.gapArray = ss.gapArray;

	if(ss.strand == 0) {
		makeMatchStringForSite(ss, read, sizes, sites, r, whole_genome, threadId);
	}
	else {
		//cout << read_reverse << endl;
		makeMatchStringForSite(ss, read_reverse, sizes, sites, r, whole_genome, threadId);
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
		int max_score, int max_quick_score, bool fully_defined, int strand, string &whole_genome) {
	vector<int> starts;
	vector<int> stops;
	int get_hits = getHits(keys, MAX_LEN, starts, stops, sizes);

	/*for(unsigned int i = 0; i < starts.size(); i++) {
		for(unsigned int j = starts[i]; j < stops[i]; j++) {
			cout << sites[j] << " ";
		}
	}
	cout << endl;*/
	/*for(unsigned int i = 0; i < stops.size(); i++) {
		cout << stops[i] << " ";
	}
	cout << endl;
*/
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

	Heap heap = Heap();

	vector<int> sites_tmp;
	heap.clear();
	heap.initHeap(offsets.size());

	for(int i = 0; i < get_hits; i++){
		int start = starts[i];

		int a = sites[start];
		int a2 = a - offsets[i];

		Triplet t(i, start, a2);
		sites_tmp.push_back(a2);

		heap.add(t);
	}
	if(write_output) cout << "\n";

	vector<int> loc_array;
	SiteScore prevSS;

	int numberOfSites = 0;



	while(!heap.isEmpty()) {
		Triplet t = heap.peek();
		//cout << "left: " << t.site << endl;
		/*for(unsigned int i = 0; i < heap.array.size(); i++) {
			cout << heap.array[i].site << endl;
		}*/
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
					score = extendScore(read, offsets, sites_tmp, center_index, loc_array, get_hits, approx_hits, whole_genome);

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
						results[results.size()-1].score = betterScore;
						results[results.size()-1].perfect = (prevSS.perfect || perfect1 || perfect2);
						//if(prevSS.perfect){prevSS.semiperfect=true;}

					} else if(shortEnough && prevSS.start == site2){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[results.size()-1].stop = site3;
							if(!prevSS.perfect) { perfects_found++; }
							results[results.size()-1].perfect = true;
							results[results.size()-1].score = betterScore;
						}else{
							results[results.size()-1].stop = maxStop;
							// prevSS.setPerfect(bases); TODO add perfect set
						}
						results[results.size()-1].score = betterScore;
					}else if(shortEnough && prevSS.stop == site3){
						if(perfect2){
							//do nothing
						}else if(perfect1){
							results[results.size()-1].start = site2;
							if(!prevSS.perfect) { perfects_found++; }
							results[results.size()-1].perfect = true;
							results[results.size()-1].score = betterScore;
						}else{
							results[results.size()-1].start = minStart;
							//prevSS.setPerfect(bases); TODO add perfect set
						}
						results[results.size()-1].score = betterScore;
					}else if(shortEnough && (maxStop - minStart <= (int) read.size() + MAX_SUBSUMPTION_LENGTH) && !perfect1 && !perfect2){
						results[results.size()-1].start = minStart;
						results[results.size()-1].stop = maxStop;
						results[results.size()-1].score = betterScore;
						//prevSS.setPerfect(bases); TODO add perfect set
					}else{
						ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
						//cout << "first: start: " << site2 << ", stop: " << site3 << "\n";
						//cout << "length: " << (site3 - site2) << endl;
						//cout << "score: " << score << endl;
						//cout << "size of results: " << results.size() << endl;
						//if(!perfect1){ss.setPerfect(bases);} TODO add perfect set
					}

				}else if(inbounds){
					ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
					//cout << "second: start: " << site2 << ", stop: " << site3 << "\n";
					//cout << "length: " << (site3 - site2) << endl;
					//cout << "score: " << score << endl;
					//cout << "size of results: " << results.size() << endl;
					//if(!perfect1){ss.setPerfect(bases);} TODO add perfect setter
				}

				if(isNull(ss) == false){
					results.push_back(ss);
					numberOfSites++;
					if(ss.perfect){
						if(isNull(prevSS) == true || !prevSS.perfect || overlap(ss.start, ss.stop, prevSS.start, prevSS.stop)){
							perfects_found++;
							if(perfects_found >= 2) {
								heap.clear();
								break;
							}
						}
					}
					prevSS = ss;
				}
			}
		}

		while(heap.peek().site == site){

			Triplet t2 = heap.poll();
			int row = t2.row+1;
			int col = t2.column;
			if(row < stops[col]){
				t2.row=row;
				int a = sites[row];
				int a2 = a-offsets[col];
				t2.site = a2;
				sites_tmp[col]=a2;
				heap.add(t2);
			}
			else if(heap.size() < approx_hits_cutoff){
				/*cout << "--------------" << endl;
				for(unsigned int i = 0; i < results.size(); i++) {
					cout << "start: " << results[i].start << ", stop: " << results[i].stop << endl;
				}*/
				return;
			}
		}
	}
}

void processRead(int *sizes, int *sites, string &r1, vector<Result> &resultsFinal, string &whole_genome, int threadId, int br) {

	//cout << " " << br << " ";

	int split_size = r1.size() / split_count;
	for(int i = 0; i < split_count; i++) {
		/*string read;
		if(i == 3) {
			read = r1.substr(1450,550);
		}
		else {
			read = r1.substr(i * 500,500);
		}*/
		string read = r1.substr(i * split_size, split_size);
	/*	string r2 = read.substr(500,500);
		string r3 = read.substr(1000,500);
		string r4 = read.substr(1500,500);*/

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

		//cout << "Read: " << br << " prepared." << endl;

		if(prescan_qscore){
			vector<vector<int> > prescanResults;
			prescanAllBlocks(prescanResults, bestScores, read_keys_final, offsets,
					keys_reversed, offsets_reversed, pretend_all_bases_covered, sizes, sites);

			//cout << "Read: " << br << "prescan all blocks. " << endl;

			precounts=prescanResults[0];
			prescores=prescanResults[1];

			if(bestScores[1]<MIN_APPROX_HITS_TO_KEEP) {
				//cout << "best scores: " << bestScores[1] << endl;
				return;
			}
			if(bestScores[3]<quick_max_score*MIN_QSCORE_MULT2) {
				//cout << "best scores3: " << bestScores[3] << endl;
				return;
			} //if(bestScores[3]<maxQuickScore(offsetsP, keyScoresP)*.10f){return result;}

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

		//cout << "Read: " << br << " prescaned." << endl;

	//	cout << "precounts[0]: " << precounts[0] << endl;
	//	cout << "prescores[0]: " << prescores[0] << endl;
	//	cout << "precounts[1]: " << precounts[1] << endl;
	//	cout << "prescores[1]: " << prescores[1] << endl;

		if(precounts[0] >= hitsCutoff && prescores[0] >= qscoreCutoff) {
			align(bestScores, read_keys_final, read, offsets, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 0, whole_genome);
		}
		if(precounts[1] >= hitsCutoff && prescores[1] >= qscoreCutoff) {
			align(bestScores, keys_reversed, read_reverse, offsets_reversed, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 1, whole_genome);
		}

		//cout << "Read: " << br << " aligned." << endl;
	/*
		for(unsigned int i = 0; i < results.size(); i++) {
			cout << "start: " << results[i].start << ", stop: " << results[i].stop << endl;
		}
	*/

		if(results.size() != 0) {
			makeMatchString(results, read, read_reverse, sizes, sites, resultsFinal, whole_genome, threadId, br);
			aligned_base_num += read.size();
		}

	}
	return;
}

void sortResults(vector<Result> &results) {
	for(unsigned int i = 0; i < results.size()-1; i++) {
		unsigned int min = i;
		for(unsigned int j = i+1; j < results.size(); j++) {
			if(results[j].br < results[min].br) min = j;
		}
		if(min != i) {
			Result r = results[i];
			results[i] = results[min];
			results[min] = r;
		}
	}
}

void writeResults(vector<vector<Result> > &set_of_results, map<int, Result> correct_results, string infile) {

	vector<Result> results;
	for(unsigned int i = 0; i < set_of_results.size(); i++) {
		for(unsigned int j = 0; j < set_of_results[i].size(); j++) {
			results.push_back(set_of_results[i][j]);
		}
	}
	cout << "size of results: " << results.size() << endl;
	sortResults(results);
	ofstream out_res(infile.c_str());
	int br = 1;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		map<int, Result>::iterator pos = correct_results.find(r.br);
		Result r2 = pos->second;
		out_res << r.br << "-" << r.start << "-" << r.stop << "-" << (r.stop-r.start) << endl;
		//out_res << r.precise_start << "-" << r.precise_stop << endl;
		out_res << r2.br << "-" << r2.start << "-" << r2.stop << "-" << (r2.stop-r2.start) << endl;
		out_res << r.matchString << endl;
		br++;
	}

	out_res.close();
}

int calculateGappedPrecisePercentage(int precise_start, int precise_stop, int start1, int stop1, int start2, int stop2, string &matchString) {
	int number = 0;
	int length = precise_stop - precise_start;

	int start = start2;
	if(start2 < start1) start = start1;
	int stop = stop2;
	if(stop2 > stop1) stop = stop1;

	start = start - precise_start;
	stop = stop - precise_start;

	if(start < 0 && stop < 0) return 0;
	if(start > length && stop > length) return 0;
	if(start < 0 && stop > 0 && stop < length) start = 0;
	if(start > 0 && start < length && stop > length) stop = length;
	if(start < 0 && stop > length) {
		start = 0;
		stop = length;
	}

	for(int j = start; j <= stop; j++) {
		if(matchString[j] == 'm') number++;
	}

	return number;
}

int calculateOverlap(int start1, int stop1, int start2, int stop2) {
	if(start1 <= start2 && stop1 >= stop2) return (stop2-start2);
	else if(start2 <= start1 && stop2 >= stop1) return (stop1-start1);
	else if(start1 <= start2 && stop2 >= stop1) return (stop1-start2);
	else if(start2 <= start1 && stop1 >= stop2) return (stop2-start1);
	else return 0;
}

void calculateStatistics2(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
	int sum = 0;
	if(r1.gapArray.size() == 0) {
		r1.gapArray.push_back(r1.start);
		r1.gapArray.push_back(r1.stop);
	}
	for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
		for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
			if(overlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1])) {
				int tmp = calculateOverlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1]);
				sum += tmp;
				//out_stat << r1.gapArray[j] << "-" << r1.gapArray[j+1] << "  " << r2.gapArray[i] << "-" << r2.gapArray[i+1] << " (" << tmp << ")" << endl;
			}
		}
	}
	//out_stat << "bases overlap: " << sum << endl;
	/*for(unsigned int i = 0;  i < r1.gapArray.size(); i++) {
		out_stat << r1.gapArray[i] << " ";
	}
	out_stat << endl;
	for(unsigned int i = 0;  i < r2.gapArray.size(); i++) {
		out_stat << r2.gapArray[i] << " ";
	}*/
//	out_stat << endl;
	statistics.covered += sum;

}

void calculateStatistics(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
/*
	if(r1.br == 26) {
		cout << "ere" << endl;
	}*/
	statistics.reads_with_result++;
	if(r2.gapArray.size() > 2) {
		for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				if(r1.matchString[j] == 'm') statistics.gapped_percentage++;
			}
		}

		for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
				if(overlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1])) {
					statistics.gapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1], r1.matchString);
				}
			}
		}

		statistics.gapped++;
		if(r1.start == r2.start) statistics.gapped_with_start++;
		if(r1.stop == r2.stop) statistics.gapped_with_stop++;
		if(r1.start == r2.start && r1.stop == r2.stop) statistics.gapped_with_start_and_stop++;
		if(r1.precise_start == r2.start) statistics.gapped_with_precise_start++;
		if(r1.precise_stop == r2.stop) statistics.gapped_with_precise_stop++;
		if(r1.precise_start == r2.start && r1.precise_stop == r2.stop) statistics.gapped_with_precise_start_and_stop++;
		if(r1.gapArray.size() == r2.gapArray.size()) {
			for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
				statistics.exons++;
				statistics.correct_exons++;
				statistics.same_numbered_exons++;
				if(r1.gapArray[i] == r2.gapArray[i]) statistics.correct_exon_start++;
				if(r1.gapArray[i+1] == r2.gapArray[i+1]) statistics.correct_exon_stop++;
				if(r1.gapArray[i] == r2.gapArray[i] && r1.gapArray[i+1] == r2.gapArray[i+1]) statistics.correct_exon_start_and_stop++;
			}
		}
		else {
			statistics.exons = statistics.exons + r1.gapArray.size()/2;
			statistics.correct_exons = statistics.correct_exons + r2.gapArray.size()/2;
		}
	}
	else {
		statistics.ungapped++;
		if(r1.start == r2.start) statistics.ungapped_with_start++;
		if(r1.stop == r2.stop) statistics.ungapped_with_stop++;
		if(r1.start == r2.start && r1.stop == r2.stop) statistics.ungapped_with_start_and_stop++;
		if(r1.precise_start == r2.start) statistics.ungapped_with_precise_start++;
		if(r1.precise_stop == r2.stop) statistics.ungapped_with_precise_stop++;
		if(r1.precise_start == r2.start && r1.precise_stop == r2.stop) statistics.ungapped_with_precise_start_and_stop++;
		for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				if(r1.matchString[j] == 'm') statistics.ungapped_percentage++;
			}
		}
		if(overlap(r2.gapArray[0], r2.gapArray[1], r1.precise_start, r1.precise_stop)) {
			statistics.ungapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[0], r2.gapArray[0+1], r1.precise_start, r1.precise_stop, r1.matchString);
		}

	}
	statistics.percentage = statistics.gapped_percentage + statistics.ungapped_percentage;
	statistics.precise_percentage = statistics.gapped_precise_percentage + statistics.ungapped_precise_percentage;

}

void writeStatistics(vector<vector<Result> > &set_of_results, map<int, Result> &correct_results, string infile, Statistic &statistic) {
	ofstream out_stat(infile.c_str());
	//Statistic statistic = Statistic();
	//cout << "number: " << results.size();
	for(unsigned int j = 0; j < set_of_results.size(); j++) {
		vector<Result> results = set_of_results[j];
		for(unsigned int i = 0; i < results.size(); i++) {
			/*if(results[i].br % 100 == 0) */ //cout << "stat: " << results[i].br << endl;
			Result r1 = results[i];
			map<int, Result>::iterator pos = correct_results.find(r1.br);
			Result r2 = pos->second;
			calculateStatistics(r1, r2, out_stat, statistic);
			calculateStatistics2(r1, r2, out_stat, statistic);
		}
	}

	double percentage = statistic.percentage / (double) statistic.reads_with_result;
	double precise_percentage = statistic.precise_percentage / (double) statistic.reads_with_result;
	double gapped_percentage = statistic.gapped_percentage / (double) statistic.gapped;
	double gapped_precise_percentage = statistic.gapped_precise_percentage / (double) statistic.gapped;
	double ungapped_percentage = statistic.ungapped_percentage / (double) statistic.ungapped;
	double ungapped_precise_percentage = statistic.ungapped_precise_percentage / (double) statistic.ungapped;

	out_stat << "reads: " << statistic.reads << endl;
	out_stat << "reads_with_result: " << statistic.reads_with_result << endl;
	out_stat << "gapped: " << statistic.gapped << endl;
	out_stat << "ungapped: " << statistic.ungapped << endl;
	out_stat << "gapped_with_start: " << statistic.gapped_with_start << endl;
	out_stat << "gapped_with_stop: " << statistic.gapped_with_stop << endl;
	out_stat << "gapped_with_start_and_stop: " << statistic.gapped_with_start_and_stop << endl;
	out_stat << "gapped_with_precise_start: " << statistic.gapped_with_precise_start << endl;
	out_stat << "gapped_with_precise_stop: " << statistic.gapped_with_precise_stop << endl;
	out_stat << "gapped_with_precise_start_and_stop: " << statistic.gapped_with_precise_start_and_stop << endl;
	out_stat << "ungapped_with_start: " << statistic.ungapped_with_start << endl;
	out_stat << "ungapped_with_stop: " << statistic.ungapped_with_stop << endl;
	out_stat << "ungapped_with_start_and_stop: " << statistic.ungapped_with_start_and_stop << endl;
	out_stat << "ungapped_with_precise_start: " << statistic.ungapped_with_precise_start << endl;
	out_stat << "ungapped_with_precise_stop: " << statistic.ungapped_with_precise_stop << endl;
	out_stat << "ungapped_with_precise_start_and_stop: " << statistic.ungapped_with_precise_start_and_stop << endl;
	out_stat << "exons: " << statistic.exons << endl;
	out_stat << "correct_exons: " << statistic.correct_exons << endl;
	out_stat << "same_numbered_exons: " << statistic.same_numbered_exons << endl;
	out_stat << "correct_exon_start: " << statistic.correct_exon_start << endl;
	out_stat << "correct_exon_stop: " << statistic.correct_exon_stop << endl;
	out_stat << "correct_exon_start_and_stop: " << statistic.correct_exon_start_and_stop << endl;
	out_stat << "correct_bases: " << statistic.percentage << endl;
	out_stat << "percentage: " << percentage << endl;
	out_stat << "precise_correct_bases: " << statistic.precise_percentage << endl;
	out_stat << "precise_percentage: " << precise_percentage << endl;
	out_stat << "ungapped_correct_bases: " << statistic.ungapped_percentage << endl;
	out_stat << "ungapped_percentage: " << ungapped_percentage << endl;
	out_stat << "ungapped_precise_correct_bases: " << statistic.ungapped_precise_percentage << endl;
	out_stat << "ungapped_precise_percentage: " << ungapped_precise_percentage << endl;
	out_stat << "gapped_correct_bases: " << statistic.gapped_percentage << endl;
	out_stat << "gapped_percentage: " << gapped_percentage << endl;
	out_stat << "gapped_precise_correct_bases: " << statistic.gapped_precise_percentage << endl;
	out_stat << "gapped_precise_percentage: " << gapped_precise_percentage << endl;
	out_stat << "covered: " << statistic.covered << endl;
	out_stat << "covered by aligned reads: " << 100 * (statistic.covered / (double) aligned_base_num) << "%" << endl;
	out_stat << "covered by total reads: " << 100 * (statistic.covered / (double) total_base_num) << "%" << endl;
	out_stat.close();

}


void *preProcessRead(void *threadid) {
	ThreadData3 *td = (ThreadData3 *) threadid;
	for(int i = td->start; i < td->stop; i++) {

		string read = (*(td->reads))[i];
		/*string r1 = read.substr(0,500);
		string r2 = read.substr(500,500);
		string r3 = read.substr(1000,500);
		string r4 = read.substr(1500,500);*/

		if((i+1)%100 == 0) cout << "Read: " << i+1 << " started." << endl;
		processRead(td->sizes, td->sites, read, *(td->results), *(td->whole_genome),  td->thread_id, i+1);
		/*processRead(td->sizes, td->sites, r1, *(td->results), *(td->whole_genome),  td->thread_id, i+1);
		processRead(td->sizes, td->sites, r1, *(td->results), *(td->whole_genome),  td->thread_id, i+1);
		processRead(td->sizes, td->sites, r1, *(td->results), *(td->whole_genome),  td->thread_id, i+1);*/
	}
	pthread_exit(NULL);
}

void checkParameter(string command, string value) {
	if(strcmp(command.c_str(), "-t") == 0) {
			thread_num = atoi(value.c_str());
			cout << "Number of threads se to " << thread_num << "." << endl;
		}
		else if(strcmp(command.c_str(), "-i") == 0) {
			index_location = value;
			cout << "Index location set to '" << index_location << "'." << endl;
		}
		else if(strcmp(command.c_str(), "-b") == 0) {
			build_number = atoi(value.c_str());
			cout << "Build number set to " << value << "." << endl;
		}
		else if(strcmp(command.c_str(), "-c") == 0) {
			overwrite_index = true;
			build_number = atoi(value.c_str());
			cout << "Index will be created. If index was previously created at ";
			cout << "set location it will be overwriten." << endl;
			cout << "Build number overwriten to " << value << "." << endl;
		}
		else if(strcmp(command.c_str(), "-s") == 0) {
			split_count = atoi(value.c_str());
			cout << "Reads will be split at " << value << " parts." << endl;
		}
		else if(strcmp(command.c_str(), "-k") == 0) {
			KEYLEN = atoi(value.c_str());
			cout << "KEYLEN set to " << value << "." << endl;
		}
		else {
			cout << "Parameter " << command << " not recognized." << endl;
			exit(-1);
		}
}

int main(int argc, char *argv[]) {

	if(argc != 4 && argc != 6 && argc != 8 && argc != 10) {
		cout << "Program must run with arguments: <destination_folder> <reads_file> <genome_reference_file> -t <thread_number> -i <index_location> -c <create_index>" << endl;
		cout << "<destination_folder> - folder where results will be stored." << endl;
		cout << "<reads_file> - file with reads in fasta format." << endl;
		cout << "<genome_reference_file> - file with regerence genome." << endl;
		//cout << "<create_index> - 0 if index was created before, 1 if you want to create an index for this run." << endl;
		//cout << "<reults_file> - file where resulting alignments will be stored." << endl;
		//cout << "<statistics_file> - file where statistics will be stored." << endl;
		cout << "-t <thread_number> - Optional parameter. Number of threads, default 4." << endl;
		cout << "-k <KEYLEN> - Optional parameter. Length of the key. Default 13." << endl;
		cout << "-b <build_number> - Optional parameter. Build number of the index. By default set to 1." << endl;
		cout << "-i <index_location> - Optional parameter. Location of created index if index exists at <index_location>. If index does not exist, program creates index at <index_location>. " << endl;
		cout << "By default program looks for index in the <destination_folder> location. If there is no index, by default program creates new index in <destination_folder>" << endl;
		cout << "If you wish to create index even if it exists use -c <create_index> option." << endl;
		cout << "-o <build_number> - Optional parameter. Use if you wish to create index even if index exist at given location. This will overwrite previously created index at given location. It should be called with build number; for example: -c 2." << endl;
		cout << "-s <split_count> - Optional parameter. At how many party every read should be split. Default no split.";
		cout << "" << endl;
		return -1;
	}
	struct stat sb;
	cout << "Program started: " << endl;

	if (stat(argv[1], &sb) == 0 && S_ISDIR(sb.st_mode)) {
		outdir = argv[1];
		cout << "Results will be stored in '" << outdir << "'." << endl;
		//exit(0);
	}
	else {
		outdir = argv[1];
		const string out = "mkdir -p " + outdir;
		//cout << out << endl;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << outdir << "'." << endl;
			exit(-1);
		}
		else {
			cout << "Directory '" << outdir << "' created." << endl;
			cout << "Results will be stored in '" << outdir << "'." << endl;
		//	exit(0);
		}
	}

	index_location = outdir;

	if(argc == 6) {
		checkParameter(argv[4], argv[5]);
	}
	if(argc == 8) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
	}
	if(argc == 10) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
		checkParameter(argv[8], argv[9]);
	}

	bool read_index = true;
	string genome_ref = argv[3];
	string whole_genome;

	if (stat(index_location.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		FILE *pFile;
		FILE *pFile2;
		FILE *pFile3;
		string build = SSTR(build_number);
		string sizes_location = index_location + "//sizes" + build;
		string sites_location = index_location + "//sites" + build;
		//cout << sites_location << endl;
		string info_location = index_location + "//index" + build + ".info";
		pFile = fopen( sizes_location.c_str() , "rb" );
		pFile2 = fopen(sites_location.c_str(), "rb");
		pFile3 = fopen(info_location.c_str(), "rb");

		if(!pFile || !pFile2 || !pFile3) {
			cout << "Index does not exist at location '" << index_location << "'." << endl;
			cout << "Index will be created at location '" << index_location << "'." << endl;
			read_index = false;
		}
		else {
			/*string line
			ifstream index_info(info_location.c_str());
			vector<string> elems;
			while(getline(index_info, line)) {
				split(line, ':', elems);
				if()
			}*/
		}
		if(pFile) {
			fclose(pFile);
		}
		if(pFile2) {
			fclose(pFile2);
		}
	}
	else {
		cout << "Directory '" << index_location << "' does not exist." << endl;
		read_index = false;
		const string out = "mkdir -p " + index_location;
		//cout << out << endl;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << index_location << "'." << endl;
			exit(-1);
		}
		else {
			cout << "Directory '" << index_location << "' created." << endl;
			cout << "Index will be stored in '" << index_location << "'." << endl;
		//	exit(0);
		}
	}

	if(overwrite_index) {
		read_index = false;
	}

	int **res;
	if(read_index) {
		cout << "Reading index..." << endl;
		res = readIndex(whole_genome, genome_ref, index_location);
	}
	else {
		cout << "Creating index..." << endl;
		res = createIndex(true, whole_genome, genome_ref, index_location, KEYLEN, keyspace, build_number);
	}

	cout << "genome size: " << whole_genome.size() << endl;

	int *sizes = res[1];
	int *sites = res[3];

	length_of_sizes = res[0][0];
	length_of_sites = res[2][0];

	vector<string> reads;
	map<int, Result> correct_results;
	readReads(reads, correct_results, argv[2]);

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	cout << "Number of reads: "<< reads.size() << "\n";

	pthread_t threads[thread_num];
	int rc;
	pthread_attr_t attr;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//vector<Result> results;
	vector<vector<Result> > set_of_results;
	for(int i = 0; i < thread_num; i++) {
		vector<Result> tmp;
		set_of_results.push_back(tmp);
	}
	ThreadData3 datas3[thread_num];
	int difference = reads.size() / thread_num;

	for(int i = 0; i < thread_num; i++) {
		int start =  i*difference;
		int stop =  i*difference + difference;
		//cout << "thread: " << i << " " << start << "-" << stop << endl;
		datas3[i] = ThreadData3(i, sizes, sites, &reads, start, stop, &whole_genome, &set_of_results[i]);
	}
	datas3[thread_num-1].stop = reads.size();
	for(int i = 0; i < thread_num; i++) {
		cout << "Thread " << i << " started." << endl;
		rc = pthread_create(&threads[i], NULL, preProcessRead, (void *) &datas3[i]);
	}

	pthread_attr_destroy(&attr);
	for(int i=0; i < thread_num; i++ ){
		rc = pthread_join(threads[i], &status);
		if (rc){
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id :" << i ;
		cout << "  exiting with status :" << status << endl;
	}
/*
	for(unsigned int i = 0; i < reads.size(); i++) {
		//cout << (i+1) << "------------------- \n";
		processRead(sizes, sites, reads[i], results, whole_genome);
	}*/
	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Reads aligned: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	writeResults(set_of_results, correct_results, outdir + "//" + "results.txt");
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Results writen: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	Statistic s = Statistic(reads.size());
	writeStatistics(set_of_results, correct_results, outdir + "//" + "statistics.txt", s);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Statistics written: " << timefinal << " seconds. " << endl;

	if(read_index) {
		free(sizes);
		free(sites);
	}
	else {
		delete [] sizes;
		delete [] sites;
	}

	delete [] res[0];
	delete [] res[2];
	delete [] res;

	cout << "Program ended." << endl;


	return 0;
}

