#include "config.h"
#include "genome_index.h"
#include "heap.h"
#include "const.h";

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
	int start30;
	int stop30;
	int start_and_stop30;
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
		start30 = 0;
		stop30 = 0;
		start_and_stop30 = 0;
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

struct ThreadData3 {
	int thread_id;
	int *sizes;
	int *sites;
	vector<Read> *reads;
	int start;
	int stop;
	string *whole_genome;
	vector<Result> *results;
	vector<Read> *unaligned_reads;
	ThreadData3() {
	}
	ThreadData3(int thread_id_, int *sizes_, int *sites_, vector<Read> *reads_, int start_, int stop_, string *whole_genome_, vector<Result> *results_, vector<Read> *unaligned_reads_) {
		thread_id = thread_id_;
		sizes = sizes_;
		sites = sites_;
		reads = reads_;
		start = start_;
		stop = stop_;
		whole_genome = whole_genome_;
		results = results_;
		unaligned_reads = unaligned_reads_;
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

/*
 *
 * READ FILES FUNCTIONS
 *
 */

void readUnalignedReads(vector <Read> &tmp_unaligned_reads, vector <Read> &reads, string infile) {
	ifstream ifs(infile.c_str());

	ofstream out("test.txt");

	string line;
	string line2;
	string line3;
	string line4;
	vector<int> brs;

	while(getline(ifs, line)) {
		getline(ifs, line2);
		getline(ifs, line3);
		getline(ifs, line4);

		char delim2 = '-';
		stringstream ss(line);
		string item;
		int br;
		while (std::getline(ss, item, delim2)) {
			br = atoi(item.c_str());
			out << br << endl;;
			break;
		}
		brs.push_back(br);
	}
	cout << "done" << endl;

	for(unsigned int i = 1; i <= reads.size(); i++) {
		if(std::find(brs.begin(), brs.end(), i) != brs.end()) {
			/* v contains x */
		} else {
			Read tm = reads[i-1];
			tmp_unaligned_reads.push_back(tm);
		}
	}
}

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

void readFastaReads(vector<Read> &reads, string infile, map<int, FastaRead> fastaReads) {
	ifstream ifs(infile.c_str());
	if(!ifs) {
		cout << "File " << infile << " does not exist." << endl;
		exit(-1);
	}
	string line;
	string line2;
	int br = 1;
	while(getline(ifs, line)) {
		getline(ifs, line2);
		Read r(br, line2);
		reads.push_back(r);
		string first = line.substr(2, line.size()-2);
		string second = line2;
		FastaRead fr(first, second);
		std::map<int,FastaRead>::iterator it = fastaReads.begin();
		fastaReads.insert(it, pair<int,FastaRead>(br, fr));
		br++;
	}

}

void readReads(vector<Read> &reads, map<int, FastaRead> &fastaReads, map<int, Result> &results, string infile) {
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
		string first = SSTR(br);
		string second = line6;
		FastaRead fr(first, second);
		std::map<int,FastaRead>::iterator it2 = fastaReads.begin();
		fastaReads.insert(it2, pair<int,FastaRead>(br, fr));
		Result r = Result(br, atoi(line2.c_str()), atoi(line3.c_str()));
		r.matchString = line4;
		gapArrayFromString(r.gapArray, line5);
		std::map<int,Result>::iterator it = results.begin();
		results.insert(it, pair<int,Result>(br, r));
		Read rea(br, line6);
		reads.push_back(rea);
		total_base_num += line6.size();
		br++;
	}
}

/*
 *
 * FILTER REF
 *
 */

void writeNewRef(string file, string file2, vector<int> &coverage, string &part_genome, string &whole_genome, vector<int> &positions) {

	ofstream out(file.c_str());
	ofstream out2(file2.c_str());
	vector<int> starts;
	vector<int> stops;
	int start = -1;
	int stop = -1;
	int position = -1;
	for(unsigned int i = 0; i < coverage.size(); i++) {
		if(coverage[i] > 0) {
			start = i;
			stop = i;
			position = i;
			break;
		}
	}
	for(unsigned int i = position; i < coverage.size(); i++) {
		if(coverage[i] > 0) {
			//sum++;
			if(stop + 32000 > i) {
				stop = i;
			}
			else {
				starts.push_back(start-100);
				stops.push_back(stop + 100);
				start = i;
				stop = i;
			}
			/*int tmp = i - stop;

			if(tmp > 32000) {
				cout << stop << "-" << i << " " << tmp << endl;
				sum += tmp;
			}
			stop = i;*/
		}
	}
	starts.push_back(start);
	stops.push_back(stop);

	for(unsigned int i = 0; i < starts.size(); i++) {
		out << starts[i] << "-" << stops[i] << /*" " << (stops[i] - starts[i]) << */endl;
		//int pos2 = starts[i];
		int pos3 = starts[i]-32000;
		int pos4 = stops[i]+1;
		//int length = (stops[i]-starts[i])+1;
		out2 << "> " << starts[i] << "-" << stops[i] << endl;
		for(int j = 0; j < 32000; j++) {
			positions.push_back(pos3);
			part_genome += "N";
			out2 << "N";
			pos3++;
		}
		for(int j =  starts[i]; j <= stops[i]; j++) {
			positions.push_back(j);
			out2 << whole_genome[j];
			part_genome.push_back(whole_genome[j]);
		}
		//out2 << whole_genome.substr(pos2, length);
		for(int j = 0; j < 32000; j++) {
			positions.push_back(pos4);
			out2 << "N";
			part_genome.push_back('N');
			pos4++;
		}
		out2 << endl;
	}

}

void writeNewRef2(string file, vector<int> &coverage) {
	int start = 0;
	int stop = 0;
	int previous = coverage[0];
	vector<int> starts;
	vector<int> stops;
	for(unsigned int i = 1; i < coverage.size(); i++) {
		if(coverage[i] != 0 && previous != 0) {
			stop = i;
			previous = coverage[i];
		}
		else if(coverage[i] != 0 && previous == 0) {
			start = i;
			stop = i;
			previous = coverage[i];
		}
		else if(coverage[i] == 0 && previous != 0) {
			starts.push_back(start);
			stops.push_back(stop);
			start = i+1;
			stop = i+1;
			previous = coverage[i];
		}
		else if(coverage[i] == 0 && previous == 0) {
			start = i+1;
			stop = i+1;
			previous = coverage[i];
		}

	}
}

/*
 *
 * COVERAGE
 *
 */


void calculateCoverageFromCig(vector<int> &coverage, string infile) {
	ifstream res(infile.c_str());
	string line;
	string linestr;
	while(getline(res, line)) {
		getline(res, line);
		getline(res, line);
		getline(res, line);
		getline(res, line);
		vector<int> array;
		gapArrayFromString(array, linestr);
		for(unsigned int k = 0; k < array.size(); k+=2) {
			int start = array[k];
			int stop = array[k+1];
			for(int l = start; l <= stop; l++) {
				coverage[l] = coverage[l]+1;
			}
		}

	getline(res, line);
	}
}


void calculateCoverageFromResults(vector<int> &coverage, string infile) {
	ifstream res(infile.c_str());
	string line;
	string linestr;
	while(getline(res, line)) {
		getline(res, line);
		getline(res, line);
		linestr = line.substr(0, line.size()-2);
		vector<int> array;
		gapArrayFromString(array, linestr);
		for(unsigned int k = 0; k < array.size(); k+=2) {
			int start = array[k];
			int stop = array[k+1];
			for(int l = start; l <= stop; l++) {
				coverage[l] = coverage[l]+1;
			}
		}

	getline(res, line);
	}
}

void calculateCoverageFromVector(vector<int> &coverage, vector<vector<Result> > &results) {
	for(unsigned int i = 0; i < results.size(); i++) {
		for(unsigned int j = 0; j < results[i].size(); j++) {
			Result r = results[i][j];
			if(r.gapArray.size() == 0) {
				r.gapArray.push_back(r.start);
				r.gapArray.push_back(r.stop);
			}
			if(r.br % 50 == 0) {
				cout << "read: " << r.br << endl;
			}
			for(unsigned int k = 0; k < r.gapArray.size(); k+=2) {
				int start = r.gapArray[k];
				int stop = r.gapArray[k+1];
				for(int l = start; l <= stop; l++) {
					coverage[l] = coverage[l]+1;
				}
			}
		}
	}
}

void writeCoverageFiltered(string path, vector <int> &coverage) {
	ofstream outs(path.c_str());
	for(unsigned int i = 0; i < coverage.size(); i++) {
	    if(coverage[i] > 0) outs << i << "\t";
	    if(coverage[i] > 0) outs << coverage[i] << endl;
	}

}

void writeCoverage(string path, vector <int> &coverage) {
	ofstream outs(path.c_str());
	for(unsigned int i = 0; i < coverage.size(); i++) {
	    outs << coverage[i] << "\t";
	}
	outs << endl;
}

void checkCoverage(vector<int> &coverage) {
	for(unsigned int i = 0; i < coverage.size(); i++) {
	    if(coverage[i] > 0) cout << coverage[i] << " ";
	}
	cout << endl;
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

	//	cout << "approxHitsCutoff: " << approxHitsCutoff << endl;
	//	cout << "approxHits: " << approxHits << endl;

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

	int bestqscore = 0;
	int maxHits = 0;
	int minHitsToScore = MIN_APPROX_HITS_TO_KEEP;

	int maxQuickScore_=calcMaxQuickScore(offsetsP, keysP);

	vector<int> counts;
	vector<int> scores;
	prescanResults.push_back(counts);
	prescanResults.push_back(scores);

	vector<int> keys = keysP;
	vector<int> offsets = offsetsP;

	vector<int> starts;
	vector<int> stops;

	int numHits = getHits(keys, std::numeric_limits<int>::max(), starts, stops, sizes);

	if(numHits < minHitsToScore){
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

		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;

		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
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

		vector<Triplet> triples;
		vector<int> values;

		vector<int> temp;
		findMaxQscore2(starts, stops, offsets, triples, values, keys, minHitsToScore, true,
				bestqscore>=maxQuickScore_ && allBasesCovered, sizes, sites, temp);

		prescanResults[0].push_back(temp[1]);
		prescanResults[1].push_back(temp[0]);

		bestqscore = max(temp[0], bestqscore);
		maxHits = max(maxHits, temp[1]);
		if(bestqscore >= maxQuickScore_ && allBasesCovered){

			minHitsToScore=max(minHitsToScore, maxHits);
			bestScores[1]=max(bestScores[1], maxHits);
			bestScores[3]=max(bestScores[3], bestqscore);
			return;
		}
	}


	bestScores[1]=max(bestScores[1], maxHits);
	bestScores[3]=max(bestScores[3], bestqscore);

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

void fillUnlimited(string &read, string &ref, int refStartLoc, int refEndLoc, long max[], vector<vector<vector<long> > > &packed) {
	int rows = read.size();
	int columns = refEndLoc-refStartLoc+1;

	long maxGain = (read.size()-1) * POINTSoff_MATCH2+POINTSoff_MATCH;
	long subfloor = 0-2*maxGain;
	long BARRIER_I2 = rows-BARRIER_I1;
	long BARRIER_I2b = columns-1;
	long BARRIER_D2 = rows-BARRIER_D1;

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

				long scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				long streak=(packed[MODE_MS][row-1][col-1]&TIMEMASK);

				{//Calculate match/sub score

					if(match){

						long scoreMS=scoreFromDiag+(prevMatch ? POINTSoff_MATCH2 : POINTSoff_MATCH);

						long scoreD=scoreFromDel+POINTSoff_MATCH;
						long scoreI=scoreFromIns+POINTSoff_MATCH;

						long score;
						long time;
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

						long scoreMS;
						if(ref1!='N' && call1!='N'){
							scoreMS=scoreFromDiag+(prevMatch ? (streak<=1 ? POINTSoff_SUBR : POINTSoff_SUB) :
								POINTSoff_SUB_ARRAY[streak+1]);
						}else{
							scoreMS=scoreFromDiag+POINTSoff_NOCALL;
						}

						long scoreD=scoreFromDel+POINTSoff_SUB; //+2 to move it as close as possible to the deletion / insertion
						long scoreI=scoreFromIns+POINTSoff_SUB;

						long score;
						long time;
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

				long streak=packed[MODE_DEL][row][col-1]&TIMEMASK;

				long scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;

				long scoreMS=scoreFromDiag+POINTSoff_DEL;

				long scoreD=scoreFromDel+(streak==0 ? POINTSoff_DEL :
					streak<LIMIT_FOR_COST_3 ? POINTSoff_DEL2 :
						streak<LIMIT_FOR_COST_4 ? POINTSoff_DEL3 :
							streak<LIMIT_FOR_COST_5 ? POINTSoff_DEL4 :
								((streak&MASK5)==0 ? POINTSoff_DEL5 : 0));

				if(ref1=='N'){
					scoreMS+=POINTSoff_DEL_REF_N;
					scoreD+=POINTSoff_DEL_REF_N;
				}else if(gap){
					long tmp = scoreMS;
					long tmp2 = scoreD;
					scoreMS+=POINTSoff_GAP;
					scoreD+=POINTSoff_GAP;
				}

				long score;
				long time;
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

				long streak=packed[MODE_INS][row-1][col]&TIMEMASK;

				long scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;

				long scoreMS=scoreFromDiag+POINTSoff_INS;
				long scoreI=scoreFromIns+POINTSoff_INS_ARRAY[streak+1];

				long score;
				long time;
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
			long x=packed[state][rows][col]&SCOREMASK;
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
/*
void traceback2(string &read, string &ref, int refStartLoc, int refEndLoc, int row, int col, int state, Result &r, vector<vector<vector<long> > > &packed) {

	int outPos=0;
	int gaps=0;

	int addon = refEndLoc - col;

	string out;
	string rout;

	while(row>0 && col>0){

		long time=packed[state][row][col]&TIMEMASK;
		long prev;

		if(state==MODE_MS){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
				else{prev=MODE_INS;}
			}

			char c=toupper(read[row-1]);
			char r=toupper(ref[refStartLoc+col-1]);
			if(c==r){
				out.push_back('M');
				rout.push_back(r);
			}else{
				out.push_back('S');
				rout.push_back(r);
			}
			row--;
			col--;
		}else if(state==MODE_DEL){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
				else{prev=MODE_DEL;}
			}

			char r=toupper(ref[refStartLoc+col-1]);
			if(r==GAPC){
				out.push_back('-');
				gaps++;
				rout.push_back('-');
			}else{
				out.push_back('D');
				rout.push_back(r);
			}
			col--;
		}else{
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
				if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else{prev=MODE_INS;}
			}

			if(col==0){
				rout.push_back('x');
			}else{
				rout.push_back('I');
				out.push_back(toupper(read[row-1]));
			}
			row--;
		}

		state=prev;
		outPos++;
	}

	if(col!=row){
		while(row>0){
			out.push_back(toupper(read[row-1]));
			rout.push_back('i');
			outPos++;
			row--;
			col--;
		}
		if(col>0){
			//do nothing
		}
	}

	reverse(out.begin(),out.end());
	reverse(rout.begin(),rout.end());
	cout << out << endl;
	cout << rout << endl;

	r.precise_start = r.start-MARIC_PADDING + col;
	r.precise_stop = r.stop+MARIC_PADDING -addon;

	if(gaps==0){
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
	return;
}
*/

void traceback(string &read, string &ref, int refStartLoc, int refEndLoc, int row, int col, int state, Result &r, vector<vector<vector<long> > > &packed) {

	long outPos=0;
	long gaps=0;

	long addon = refEndLoc - col;

	string out;

	while(row>0 && col>0){

		long time=packed[state][row][col]&TIMEMASK;
		long prev;

		if(state==MODE_MS){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row-1][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row-1][col-1]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel && scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else if(scoreFromDel>=scoreFromIns){prev=MODE_DEL;}
				else{prev=MODE_INS;}
			}

			char c=toupper(read[row-1]);
			char r=toupper(ref[refStartLoc+col-1]);
			if(c==r){
				out.push_back('m');
			}else{
				out.push_back('s');
			}
			row--;
			col--;
		}else if(state==MODE_DEL){
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row][col-1]&SCOREMASK;
				long scoreFromDel=packed[MODE_DEL][row][col-1]&SCOREMASK;
				if(scoreFromDiag>=scoreFromDel){prev=MODE_MS;}
				else{prev=MODE_DEL;}
			}

			char r=toupper(ref[refStartLoc+col-1]);
			if(r==GAPC){
				out.push_back('-');
				gaps++;
			}else{
				out.push_back('d');
			}
			col--;
		}else{
			if(time>1){prev=(long)state;}
			else{
				long scoreFromDiag=packed[MODE_MS][row-1][col]&SCOREMASK;
				long scoreFromIns=packed[MODE_INS][row-1][col]&SCOREMASK;
				if(scoreFromDiag>=scoreFromIns){prev=MODE_MS;}
				else{prev=MODE_INS;}
			}

			if(col==0){
				out.push_back('x');
			}else{
				out.push_back('i');
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

	r.precise_start = r.start-MARIC_PADDING + col;
	r.precise_stop = r.stop+MARIC_PADDING -addon;

	if(gaps==0){
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
	return;
}


void fillLimited(string &read, string &ref, int refStartLoc, int refEndLoc, int score, vector<int> &gaps, long max[], Result &r) {
	string gref;
	int grefLimit = (refEndLoc - refStartLoc) + 2 * MARIC_PADDING;
	//bool gapped = true;
	if((int)gaps.size() == 0) {
		//gapped = false;
		for(int i = refStartLoc-MARIC_PADDING; i < refEndLoc + MARIC_PADDING; i++) {
			gref.push_back(ref[i]);
		}
	}
	else {
		grefLimit = makeGref(ref, gaps, (refStartLoc - MARIC_PADDING), (refEndLoc + MARIC_PADDING), gref);
	}

	int rows = read.size();
	int columns = grefLimit+1;
	vector<vector<vector<long> > > packed;

	for(int matrix=0; matrix<3; matrix++){
		vector<long> row;
		for(int i = 0; i <= columns; i++) {
			row.push_back(0);
		}
		vector<vector<long> > mat;
		mat.push_back(row);
		for(int i=1; i<=rows; i++){
			vector<long> row;
			for(int j=0; j<columns+1; j++){
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

	fillUnlimited(read, gref, 0, grefLimit, max, packed);
	traceback(read, gref, 0, grefLimit, max[0], max[1], max[2], r, packed);
	if(r.gapArray.size() == 0) r.precise_stop--;

}

void makeMatchStringForSite(SiteScore ss, string &read, int *sizes, int *sites, Result &r, string &whole_genome, int threadId) {

	if(ss.perfect) {
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

	long max[4];
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
	}
	Result r = Result(br, ss.start, ss.stop);
	r.gapArray = ss.gapArray;
	/*if(drugi) {
		if(ss.strand == 0) {
			makeMatchStringForSite(ss, read, sizes, sites, r, whole_genome, threadId);
		}
		else {
			makeMatchStringForSite(ss, read_reverse, sizes, sites, r, whole_genome, threadId);
		}
	}*/
/*
	if(ss.strand == 0) {
		makeMatchStringForSite(ss, read, sizes, sites, r, whole_genome, threadId);
	}
	else {
		makeMatchStringForSite(ss, read_reverse, sizes, sites, r, whole_genome, threadId);
	}*/

	resultsFinal.push_back(r);
}


/*
 *
 * ALIGNER
 *
 */

void printLocationArray(vector<int> locArray) {
	for(unsigned int i = 0; i < locArray.size(); i++) {
		cout << locArray[i] << " ";
	}
	cout << endl;
}

void align(int bestScores[], vector<int> &keys, string &read, vector<int> &offsets, int *sizes,
		int *sites, vector<SiteScore> &results, bool all_bases_covered,
		int max_score, int max_quick_score, bool fully_defined, int strand, string &whole_genome) {
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

	vector<int> loc_array;
	SiteScore prevSS;

	int numberOfSites = 0;

	while(!heap.isEmpty()) {
		Triplet t = heap.peek();
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

		//cout << "approx hits: " << approx_hits << endl;
		//cout << "approx_hits_cutoff: " << approx_hits_cutoff << endl;

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

				vector<int> gapArray;

				//printLocationArray(loc_array);
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
						//printLocationArray(loc_array);
						//if(!perfect1){ss.setPerfect(bases);} TODO add perfect set
					}

				}else if(inbounds){
					ss = SiteScore(site2, site3, score, approx_hits, perfect1, strand, gapArray);
					//printLocationArray(loc_array);
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
				return;
			}
		}
	}
}

void processRead(int *sizes, int *sites, string &r1, vector<Result> &resultsFinal, vector<Read> &unaligned_reads, string &whole_genome, int threadId, int br) {

	//cout << "Alignin read: " << br << endl;

	int split_size = r1.size() / split_count;
	for(int i = 0; i < split_count; i++) {
		string read = r1.substr(i * split_size, split_size);
		string read_reverse;
		reverseComplementRead(read_reverse, read);

		vector<int> offsets;
		makeOffsets(read, offsets);
		//cout << "offsets size: " << offsets.size() << endl;
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
			vector<vector<int> > prescanResults;
			prescanAllBlocks(prescanResults, bestScores, read_keys_final, offsets,
					keys_reversed, offsets_reversed, pretend_all_bases_covered, sizes, sites);

			precounts=prescanResults[0];
			prescores=prescanResults[1];

			if(bestScores[1]<MIN_APPROX_HITS_TO_KEEP) {
				Read rea(br, read);
				unaligned_reads.push_back(rea);
				//cout << "Failed prescan 1." << endl;
				return;
			}
			if(bestScores[3]<quick_max_score*MIN_QSCORE_MULT2) {
				Read rea(br, read);
				unaligned_reads.push_back(rea);
				//cout << "Failed prescan 2." << endl;
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

		if(precounts[0] >= hitsCutoff && prescores[0] >= qscoreCutoff) {
			align(bestScores, read_keys_final, read, offsets, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 0, whole_genome);
		}
		if(precounts[1] >= hitsCutoff && prescores[1] >= qscoreCutoff) {
			align(bestScores, keys_reversed, read_reverse, offsets_reversed, sizes, sites, results, all_bases_covered, max_score, quick_max_score, fully_defined, 1, whole_genome);
		}

		if(results.size() != 0) {
			//cout << "aligned: " << br << endl;
			makeMatchString(results, read, read_reverse, sizes, sites, resultsFinal, whole_genome, threadId, br);
			aligned_base_num += read.size();
		}
		else {
			Read rea(br, read);
			unaligned_reads.push_back(rea);
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

/*
 *
 * WRITE FILES FUNCTIONS
 *
 */

string compressMatchString(string &matchString) {
	char current = matchString[0];
	int number = 1;
	string compressed;
	for(unsigned int i = 1; i < matchString.size(); i++) {
		if(current == matchString[i]) {
			number++;
		}
		else if((current == 'm' && matchString[i] == 's') || (current == 's' && matchString[i] == 'm')) {
			number++;
		}
		else {
			if(current == 's') current = 'm';
			compressed += SSTR(number);
			compressed.push_back(current);
			number = 1;
			current = matchString[i];
		}
	}
	if(current == 's') current = 'm';
	compressed += SSTR(number);
	compressed.push_back(current);
	return compressed;
}

void writeSamResults(vector<Result> &set_of_results, vector<Read> &reads, map<int, FastaRead> &read_names, string infile) {
	vector<Result> results;
	for(unsigned int j = 0; j < set_of_results.size(); j++) {
		results.push_back(set_of_results[j]);
	}

	cout << "size of results: " << results.size() << endl;
	sortResults(results);
	ofstream out_res(infile.c_str());
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		if(r.gapArray.size() == 0) {
			r.gapArray.push_back(r.start);
			r.gapArray.push_back(r.stop);
		}
		map<int, FastaRead>::iterator pos = read_names.find(r.br);
		FastaRead fr = pos->second;

		out_res << (fr.name) << "\t";
		for(unsigned j = 0; j < r.gapArray.size(); j+=2) {
			out_res << r.gapArray[j] << "-" << r.gapArray[j+1];
			if(j+2 != r.gapArray.size()) out_res << ",";
		}
		out_res << "\t";
		out_res << "0\t";
		out_res << "chr1\t";
		out_res << r.start << "\t";
		out_res << "*\t";
		out_res << compressMatchString(r.matchString) << "\t";
		out_res << "*\t";
		out_res << "0\t";
		out_res << "0\t";
		out_res << (fr.read) << "\n";
		out_res << "*\t";
		out_res << "*\t";
		out_res << "\n";
	}
}

void writeResults(vector<Result> &set_of_results, map<int, Result> &correct_results, string infile, bool second, string &whole_genome) {

	vector<Result> results;
	for(unsigned int j = 0; j < set_of_results.size(); j++) {
		results.push_back(set_of_results[j]);
	}
	cout << "size of results: " << results.size() << endl;
	sortResults(results);
	ofstream out_res(infile.c_str());
	int br = 1;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		map<int, Result>::iterator pos = correct_results.find(r.br);
		Result r2 = pos->second;
		int start = r.start;
		int stop = r.stop;
		if(second) {
			start = positions_index[start];
			stop = positions_index[stop];
			out_res << r.br << "-" << start << "-" << stop << "-" << (stop-start) << endl;
		}
		out_res << r.br << "-" << r.start << "-" << r.stop << "-" << (r.stop - r.start) << endl;
		out_res << r2.br << "-" << r2.start << "-" << r2.stop << "-" << (r2.stop-r2.start) << endl;
		if(r.gapArray.size() == 0) {
			r.gapArray.push_back(start);
			r.gapArray.push_back(stop);
		}
		for(unsigned int j = 0; j < r.gapArray.size(); j+=2) {
			out_res << r.gapArray[j] << "-" << r.gapArray[j+1] << ", ";
		}
		out_res << endl;
		//out_res << r.matchString << endl;
		out_res << whole_genome.substr(r.start, (r.stop - r.start)) << endl;
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

void calculateStatistics3(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
	if(abs(r1.start-r2.start) < 30) {
		statistics.start30++;
	}
	if(abs(r1.stop-r2.stop) < 30) {
		statistics.stop30++;
	}
	if(abs(r1.start-r2.start) < 30 && abs(r1.stop-r2.stop) < 30) {
		statistics.start_and_stop30++;
	}
}

void calculateStatistics5(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
	if(abs(positions_index[r1.start]-r2.start) < 30) {
		statistics.start30++;
	}
	if(abs(positions_index[r1.stop]-r2.stop) < 30) {
		statistics.stop30++;
	}
	if(abs(positions_index[r1.start]-r2.start) < 30 && abs(positions_index[r1.stop]-r2.stop) < 30) {
		statistics.start_and_stop30++;
	}
}

void calculateStatistics4(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
	int sum = 0;
	//cout << "Calc stat 2" << endl;

	if(r1.gapArray.size() == 0) {
		r1.gapArray.push_back(r1.start);
		r1.gapArray.push_back(r1.stop);
	}
	for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
		for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
			if(overlap(r2.gapArray[i], r2.gapArray[i+1], positions_index[r1.gapArray[j]], positions_index[r1.gapArray[j+1]])) {
				int tmp = calculateOverlap(r2.gapArray[i], r2.gapArray[i+1], positions_index[r1.gapArray[j]], positions_index[r1.gapArray[j+1]]);
				sum += tmp;
			}
		}
	}
	statistics.covered += sum;

}


void calculateStatistics2(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {
	int sum = 0;
	//cout << "Calc stat 2" << endl;

	if(r1.gapArray.size() == 0) {
		r1.gapArray.push_back(r1.start);
		r1.gapArray.push_back(r1.stop);
	}
	for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
		for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
			if(overlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1])) {
				int tmp = calculateOverlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1]);
				sum += tmp;
			}
		}
	}
	statistics.covered += sum;

}

void calculateStatistics6(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {

	//cout << "Calc stat 1" << endl,
	statistics.reads_with_result++;
	if(r2.gapArray.size() > 2) {
		/*for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			//cout << "greska -" << start << endl;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				//if(r1.matchString[j] == 'm') statistics.gapped_percentage++;
			}
		}*/

	/*	for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
				if(overlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1])) {
					statistics.gapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1], r1.matchString);
				}
			}
		}*/

		statistics.gapped++;
		if(positions_index[r1.start] == r2.start) statistics.gapped_with_start++;
		if(positions_index[r1.stop] == r2.stop) statistics.gapped_with_stop++;
		if(positions_index[r1.start] == r2.start && positions_index[r1.stop] == r2.stop) statistics.gapped_with_start_and_stop++;
	//	if(positions_index[r1.precise_start] == r2.start) statistics.gapped_with_precise_start++;
	//	if(positions_index[r1.precise_stop] == r2.stop) statistics.gapped_with_precise_stop++;
	//	if(positions_index[r1.precise_start] == r2.start && positions_index[r1.precise_stop] == r2.stop) statistics.gapped_with_precise_start_and_stop++;
	/*	if(r1.gapArray.size() == r2.gapArray.size()) {
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
		}*/
	}
	else {
		statistics.ungapped++;
		if(positions_index[r1.start] == r2.start) statistics.ungapped_with_start++;
		if(positions_index[r1.stop] == r2.stop) statistics.ungapped_with_stop++;
		if(positions_index[r1.start] == r2.start && positions_index[r1.stop] == r2.stop) statistics.ungapped_with_start_and_stop++;
	//	if(positions_index[r1.precise_start] == r2.start) statistics.ungapped_with_precise_start++;
	//	if(positions_index[r1.precise_stop] == r2.stop) statistics.ungapped_with_precise_stop++;
	//	if(positions_index[r1.precise_start] == r2.start && positions_index[r1.precise_stop] == r2.stop) statistics.ungapped_with_precise_start_and_stop++;
		/*for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			//cout << "greska 2-" << start << endl ;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				//if(r1.matchString[j] == 'm') statistics.ungapped_percentage++;
			}
		}*/
		/*if(overlap(r2.gapArray[0], r2.gapArray[1], r1.precise_start, r1.precise_stop)) {
			statistics.ungapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[0], r2.gapArray[0+1], r1.precise_start, r1.precise_stop, r1.matchString);
		}*/

	}
	statistics.percentage = statistics.gapped_percentage + statistics.ungapped_percentage;
	statistics.precise_percentage = statistics.gapped_precise_percentage + statistics.ungapped_precise_percentage;

}

void calculateStatistics(Result &r1, Result &r2, ofstream &out_stat, Statistic &statistics) {

	//cout << "Calc stat 1" << endl,
	statistics.reads_with_result++;
	if(r2.gapArray.size() > 2) {
		/*for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			//cout << "greska -" << start << endl;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				//if(r1.matchString[j] == 'm') statistics.gapped_percentage++;
			}
		}*/

	/*	for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			for(unsigned int j = 0; j < r1.gapArray.size(); j+=2) {
				if(overlap(r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1])) {
					statistics.gapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[i], r2.gapArray[i+1], r1.gapArray[j], r1.gapArray[j+1], r1.matchString);
				}
			}
		}*/

		statistics.gapped++;
		if(r1.start == r2.start) statistics.gapped_with_start++;
		if(r1.stop == r2.stop) statistics.gapped_with_stop++;
		if(r1.start == r2.start && r1.stop == r2.stop) statistics.gapped_with_start_and_stop++;
		if(r1.precise_start == r2.start) statistics.gapped_with_precise_start++;
		if(r1.precise_stop == r2.stop) statistics.gapped_with_precise_stop++;
		if(r1.precise_start == r2.start && r1.precise_stop == r2.stop) statistics.gapped_with_precise_start_and_stop++;
	/*	if(r1.gapArray.size() == r2.gapArray.size()) {
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
		}*/
	}
	else {
		statistics.ungapped++;
		if(r1.start == r2.start) statistics.ungapped_with_start++;
		if(r1.stop == r2.stop) statistics.ungapped_with_stop++;
		if(r1.start == r2.start && r1.stop == r2.stop) statistics.ungapped_with_start_and_stop++;
		if(r1.precise_start == r2.start) statistics.ungapped_with_precise_start++;
		if(r1.precise_stop == r2.stop) statistics.ungapped_with_precise_stop++;
		if(r1.precise_start == r2.start && r1.precise_stop == r2.stop) statistics.ungapped_with_precise_start_and_stop++;
		/*for(unsigned int i = 0; i < r2.gapArray.size(); i+=2) {
			unsigned int start = r2.gapArray[i] - r2.start;
			unsigned int stop = r2.gapArray[i+1] - r2.start;
			if(start > r1.matchString.size()) break;
			//cout << "greska 2-" << start << endl ;
			if(stop > r1.matchString.size()) stop = r1.matchString.size();
			for(unsigned int j = start; j <= stop; j++) {
				//if(r1.matchString[j] == 'm') statistics.ungapped_percentage++;
			}
		}*/
		/*if(overlap(r2.gapArray[0], r2.gapArray[1], r1.precise_start, r1.precise_stop)) {
			statistics.ungapped_precise_percentage += calculateGappedPrecisePercentage(r1.precise_start, r1.precise_stop, r2.gapArray[0], r2.gapArray[0+1], r1.precise_start, r1.precise_stop, r1.matchString);
		}*/

	}
	statistics.percentage = statistics.gapped_percentage + statistics.ungapped_percentage;
	statistics.precise_percentage = statistics.gapped_precise_percentage + statistics.ungapped_precise_percentage;

}

void writeStatistics(vector<Result> &results, map<int, Result> &correct_results, string infile, Statistic &statistic) {
	ofstream out_stat(infile.c_str());

	for(unsigned int i = 0; i < results.size(); i++) {
		Result r1 = results[i];
		map<int, Result>::iterator pos = correct_results.find(r1.br);
		Result r2 = pos->second;
		if(drugi == false) {
		calculateStatistics(r1, r2, out_stat, statistic);
		}
		//cout << "1" << endl;
		if(drugi) {
			calculateStatistics6(r1, r2, out_stat, statistic);
			calculateStatistics4(r1, r2, out_stat, statistic);
			calculateStatistics5(r1, r2, out_stat, statistic);
		} else {
		calculateStatistics2(r1, r2, out_stat, statistic);
		//cout << "2" << endl;
		calculateStatistics3(r1, r2, out_stat, statistic);
		//cout << "3" << endl;
		}
	}

	//cout << "gotovo" << endl;

	double percentage = statistic.percentage / (double) statistic.reads_with_result;
	double precise_percentage = statistic.precise_percentage / (double) statistic.reads_with_result;
	double gapped_percentage = statistic.gapped_percentage / (double) statistic.gapped;
	double gapped_precise_percentage = statistic.gapped_precise_percentage / (double) statistic.gapped;
	double ungapped_percentage = statistic.ungapped_percentage / (double) statistic.ungapped;
	double ungapped_precise_percentage = statistic.ungapped_precise_percentage / (double) statistic.ungapped;

	out_stat << "reads: " << statistic.reads << endl;
	out_stat << "reads_with_result: " << statistic.reads_with_result / split_count << endl;
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
	out_stat << "start30: " << statistic.start30 << endl;
	out_stat << "stop30: " << statistic.stop30 << endl;
	out_stat << "start_and_stop30: " << statistic.start_and_stop30 << endl;
	out_stat.close();

}


void *preProcessRead(void *threadid) {
	ThreadData3 *td = (ThreadData3 *) threadid;
	for(int i = td->start; i < td->stop; i++) {

		Read read = (*(td->reads))[i];

		if(drugi) {
			if((i+1)%10 == 0) cout << "Read: " << read.br << " started." << endl;
		}
		/*if((i+1)%10 == 0)*/ //cout << "Read: " << read.br << " started." << endl;
		processRead(td->sizes, td->sites, read.content, *(td->results), *(td->unaligned_reads), *(td->whole_genome),  td->thread_id, read.br);
	}
	pthread_exit(NULL);
}

/*
 *
 * PARAMS FUNCTIONS
 *
 */

void checkParameter(string command, string value) {
	if(strcmp(command.c_str(), "-t") == 0) {
		thread_num = atoi(value.c_str());
		cout << "Number of threads set to " << thread_num << "." << endl;
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
	else if(strcmp(command.c_str(), "-p") == 0) {
		precise = false;
	}
	else {
		cout << "Parameter " << command << " not recognized." << endl;
		exit(-1);
	}
}

int main3(int argc, char *argv[]) {
	vector<int> coverage;
	for(unsigned int i = 0; i < 83083540; i++) {
		coverage.push_back(0);
	//	if(i%10000 == 0) cout << i << endl;
	}
	outdir = argv[1];
	string build = "1";
	calculateCoverageFromResults(coverage, "/home/josip/Desktop/big_tests/results3fixed.txt");
	writeCoverageFiltered(outdir + "//" + "coverage3" + build + ".txt", coverage);
	//writeNewRef(outdir + "//" + "newref" + build + ".fa", coverage);
}

int main2(int argc, char *argv[]) {

	string read = "GTTAAAGCCCCAGACGTCAAAGCCATCCCGGTTGTGAATCCTTGTTTGAGAGCCTCTCCCTTCGCTCATGCTGTGTGGATCTTGAGACACCTCATCTCCCAGATGATTGTTTTGAACCTATATATTTTACCATTCTCAATGATCAAATGTAGCGCTTGAGATGGACTTACTGCTCGCTCTGCCGCATGCTGGGTTTCTTTCTAGACCCCTCTTCCTCCTGTCACCTGGGACAGCATGTGGGTTGGCAGAGGCCGATAAAGAAGCCCCGAAACTCAAATGAAAGGTCATCAGGTCCAGAGCTACCTATTAGCTTAAATTTGATCCGGTAGGAGACAAGGTGCTTAAGAGGGATAGACCTATGGGCGGAGGAAACAGCCGGAAGAAGAGAATTAAGAGACCCTCGATGAAATTTTAAATATTTAAAAAAATAAAAGTTACTGGGAAGGTGGTTAATAGAAACAGACGATTAGCGACCCAACTAGAAAGCTACACAGGTGTTCGACGAGAGAAATTGAGGGGCTGTAAGCCTAAGGATTGAATTCTAAAGATACAATGGAAAAGGGAGTATGAAGGCACACCAGGGCCTTGCGGACATCTCCTTACCTTCATTCCAGCATGAAGGAGGAGAAGGGCAGGGAGGTAGGGAGCATTTCAGCAGAGCCAGCACCAGCTTTGATAATCTCCCCCAGCTACCCTATTTAAGGAGTTCGAGGTTTAAGAGTTTAAAAACAGGTGGCACCGAGACCCCAATTCAGGAGACAGTAGGTTACTCCGGGTTCGAAGAGAAATCCTATCTGAGGCCTGTTGGCTTCCCGGGCTTCAGCAGAGACCCTCTAGGGAACCAGTAGTCACTTGATCAGTCCTGATGCCTTGACCCAATCGCCGAGAGTCGCTAAGGCTAGAGGTATTTGAGGGAGTCGAGCCCCCAAGCCTCGTGTGGCACTAATAGGCAGTTATAGAAAAGGCGGGTCCTAAGATCATTGGTGTTAGGCGGGATCGGGTAGGGGGAGCATGTCTTCTGTGCTGGAAACAGAAGGGGTATTTCAAGATGGCAGACCCATTCAATTATTGGAGCTATAAGCTCCTAGAATTGCTCGATGGGCTATCTCGGTTTCCCTTGCATCACATCTGCGCCTGAACTGCACCTGTCATGGCGGGTCCATCTTCGTCCCCGATCTCCCGTTAGTCAATCGATGTCCGCTGAACAAAATTCATTGGTGCCCCAATCACGTTCCGGTCAATCCCGCGTCTCTTCTTCCTCGGCTGGGGTTTTTCTACCACACAGCAGCGTAGTTTATGCCGAAATTCACTCTTGCTCCCTGATGGTTCCAGCTTCACGCAACTGTCTGAAGGCCCTCAGGGCTCTTTCGCTTGTGTCTGTACAACCCCCTCCACGCTCTCCCCAGATAGACCACCACCCTGGGTTCACTCATCTTCATGGGTCGTGACGGCGTGGAATGCCTGGAACCATTAACTGAGCAGTAAAGGTCAGTGGCACATAGCCCAAGAGATGAACAGGACCAGAGACAGAGGGGGGCAGGAACTCACAAGGTTATGTCTTCCTAAGAATCGGAACCCTGTCCGTGAATTCCTTCCTGAAACCTGTGTTCATTCGGTCTAGCTGTGTCCTATATAACTGAACTCTTGTGAGTACACCCCTGTGCGCACATCCGTGGCATGGCCACCTCTTCACGTATCTCCAGTACAGATGATCTATGTACAAGTTATACGTACTTGTGACAGATTTCTTCTGTCCCATGAGCAACGTTTTCCCCATCGTCCACATGTCCTATCAGTCAATGTGGTCTGACGCCCACCCTATACAGCCCAGGCGTCCCTCGCGTTGTGACCTGTGCTGCGCTCCCGCCTCACCAATTCTATGCACAAGCTTCATCAGTGAATCAGACGTAAGGAAATCTCCTAGAGTACTAGTGTCCAAGAAATAACTTCCTCCCTTCTCACTACTGCCAGCGAGGGTGGGACACATCCGGGTTTCATCT";
	string gref = "aagaacacttcagcccaaaatgtcaaaggcatctcggttgtgaaaccttgTTCTAGAGCCTCTCCCTGCTCCCGTGCTGTGAGGATCTTAAGACTCCTAAACTCACAGATGATTGTTTTCAACCTATTTATCTCACCATCCTCAATGATCAGAAGTCCACCTTGAGATGGACTTACTCCTAGTTCTGCCGCATTCTGGGTTTCTTTCCAGACCCCTCTTCCTCCTCTCACCTGGGCCAGCAGGGAGGTTGGCAGAGGGCAGGAAAGATGCCCCAAAACTTAAATGGAATGTCATCAGGTGAAGAGCTACCTCTTACCTTAACTCTTATTCGGTGGGAAACAAGGTTCCTAATAGGGATAGAGTTAAGGGCGGAGGCAACAGCAGGATGAGGATAATTAAGAGAACCTGaattaaatttttaatattttaaaaaaGGAAAGTTGCTGGGAAGCTGGTAAACAGAAACAGAAGATTAGAGGCCTAACTAGAAAGAGACACTGGTGTTAGCAGAGAGATATTGAGGGGCTGTAAGCCTAAGGTTTGAAATCTAAAGATAGAGTGGAAAAGGGAGTAGGCACCCCCACCAGCCCCTGCTTGACATCTGCTTTAGTTCATTCCAGCAGGAAGGAGGAGAAGGGCAGGGAGGTAGGGAGCTTCTCAGCAAAGCCAGCTTCAGCTTTGATAATCTCACCCACCTACCCCATTTAAGGAGTTCCAGGTTTAAGAGTTTAAAAACAGGTGGCACCCAGACCATCATTCAGGAGACAGGAACTCATTCCAGGTTCCTAGAGAACTCCTATCTCAGACCTGAGGGTTTCCAGGGCTTCAGCTGAGCTCCTCTGGCTAACCAGTAGTCACTTGATCAGTCCTGCTGCCTTGACCCCATCTCCAGGAGGGGCTATGGCCAGAGGGAGTAGAGGGAGTCCAGCCCCCAAGCCTTGTGAGGCACTGTTAGGCAGATAGGGAAAAGAGGGGTCCTTAGATCACTGGTTCAAGGAGGGATCTGGTAGGGGCAGCATTTCTTCTGGGCTGGAAACAGAATGGGGGTTTCAAGATGGCAGAACCTGAAATAAAAAAAAAGTTTTTCTTAAAGAAGATTAGCCTCAAAGAAAACCAAGGCTTTAGGAAGAAGGCACTACCTAGTGCAGACTTTAAGCTATTTCCACAGTGTGTTCTTTA------GGGTTTCTCTAATCACTCACCATTCCATTATTGGAGCTATAAGCCCCTAGAATTGCTCCATGGCCTATCTCGGTTTCCCTTGGATCTCATCTGCTCCTGAACTGCACCTGTCTGTAAAAAAACAGATGCGAGACACCTTCGTAAGTCTTCATATCCTACAGTAAGAACTCTACTTTGTGCCTTCCAGGGAAGAGGCTGACGCCATCTCCTTGACAATAGCAGCCATCTGCTAACCACCATTCCTCCTAGGGGAGTCTCAATGCCTCTTTTCACACTGGCCT----TTTCCTTGTTACCATGGCAAGTCCATCTCCGGCCCCCATCTCCCCTGAGCCAATGTGAGTCAGGTGAACAAAATTCATTGGTTCCCCAATCATGGTCCGGTCAATCCGTCTTCTCTTCTTCTTCTGCTTGGAGAAGATAAGGAGTCATGTCTTAAGGCCCCTCCTCCAGTCAGCACGGAGGAAGGGAAAAGAGTCTTCAAAACCACTGAAAAGCTAAGGGCCTGGTGGGAAGAGGGGATGTGGCATTTTAACAGGGAAATCACTAGGGCATGAAAGGCAAA---AGACTCACCGGCTGGGGTTTCTCTACCACACAGCAGCCCAGTTTGTGCCAAAATTCACTCATGTTCCCTGATGGTTCCAGCTTCACTCCGCTGTCTGAGGCCCCAAAGGGCTCTTTCCCTGGTGTTTGGACAACCCACTCCTCACTCTCCCCAGATACACCACCACCCTGGGTTCACTCAGCTGGATGGGTCCAGACAAAGTGGAATCCCTGGAACCTTTAACTGAGCAGTGAAGGTCAGTGTCTCAGAGCCTGAGAGATGAACAGGACCAGAGAGAGAGGTGGGCAGGCAGGCACAAGGTTATGTCTTCCTCAGACTCGGAACCCTGGATTAGAAGAAAGTAAAAAGAAAGATCAGTGAGAAATCTCCCCCTCCCCCAATCCTTGAGAAACCTCCTGATGAGCCTATATCTTCCCAGTCCTCCCTGAGAGGCCCTAGACCTGTACACCAACACTCATTCCCAAACACTATGAACCTCAGTTTGACTCTAAAAAAAAGGAACAGAACCAAAAGGACACCAAGGGGTTAAGGAGCAGTGGCTGGTGGTACCAAGGGAGAGGAAGGGGGCTAGAAGCACTTAGAGTTAAGAGGATGGGCGAGGACCAA----CTGTGTACAGGTTGTGGGTACATGTGCCAGGTATCTTCTGTCCCATGAGGAATATTTTCCCCACCGTCCACATGTCCTATCTGCCACTCTGCTCTGCCTCACACCCTATAGGGCCCAGGCTCCCCTCCAGCTGTGACCTGTTCTCCCCTCCCTCCTCACCAATTCTATCCTCAAGCTTCATCAGTGAACCAGACATAAGGAACTTTCCAATAGTACTAGTGTCCAAGGAAGAACTTCCTCCCTTCTCACTACTGCCAGCCAGGGTAGGACACATCTGGGTTTCTTCTTAGACT";
	int grefLimit = 2657;
	long max[4];

	int rows = read.size();
	int columns = grefLimit+1;
	vector<vector<vector<long> > > packed;

	for(int matrix=0; matrix<3; matrix++){
		vector<long> row;
		for(int i = 0; i <= columns; i++) {
			row.push_back(0);
		}
		vector<vector<long> > mat;
		mat.push_back(row);
		for(int i=1; i<=rows; i++){
			vector<long> row;
			for(int j=0; j<columns+1; j++){
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

	fillUnlimited(read, gref, 0, grefLimit, max, packed);
	Result r(0,0,0);

	traceback(read, gref, 0, grefLimit, max[0], max[1], max[2], r, packed);
}

int readParams(int argc, char *argv[]) {
	if(argc != 4 && argc != 6 && argc != 8 && argc != 10 && argc != 12 && argc != 14 && argc != 16 && argc != 18) {
		cout << "Program must run with arguments: <destination_folder> <reads_file> <genome_reference_file> -t <thread_number> -i <index_location> -c <create_index>" << endl;
		cout << "<destination_folder> - folder where results will be stored." << endl;
		cout << "<reads_file> - file with reads in fasta format." << endl;
		cout << "<genome_reference_file> - file with regerence genome." << endl;
		cout << "-t <thread_number> - Optional parameter. Number of threads, default 4." << endl;
		cout << "-k <KEYLEN> - Optional parameter. Length of the key. Default 13." << endl;
		cout << "-b <build_number> - Optional parameter. Build number of the index. By default set to 1." << endl;
		cout << "-i <index_location> - Optional parameter. Location of created index if index exists at <index_location>. If index does not exist, program creates index at <index_location>. " << endl;
		cout << "By default program looks for index in the <destination_folder> location. If there is no index, by default program creates new index in <destination_folder>" << endl;
		cout << "If you wish to create index even if it exists use -c <create_index> option." << endl;
		cout << "-c <build_number> - Optional parameter. Use if you wish to create index even if index exist at given location. This will overwrite previously created index at given location. It should be called with build number; for example: -c 2." << endl;
		cout << "-s <split_count> - Optional parameter. At how many party every read should be split. Default no split." << endl;
		cout << "-p <precision> - Optimal parameter. By default set to more precise, if set reads are aligned less precise." << endl;
		cout << "" << endl;
		return -1;
	}
	struct stat sb;
	cout << "Program started: " << endl;

	if (stat(argv[1], &sb) == 0 && S_ISDIR(sb.st_mode)) {
		outdir = argv[1];
		cout << "Results will be stored in '" << outdir << "'." << endl;
	}
	else {
		outdir = argv[1];
		const string out = "mkdir -p " + outdir;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << outdir << "'." << endl;
			return -1;
		}
		else {
			cout << "Directory '" << outdir << "' created." << endl;
			cout << "Results will be stored in '" << outdir << "'." << endl;
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
	if(argc == 12) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
		checkParameter(argv[8], argv[9]);
		checkParameter(argv[10], argv[11]);
	}
	if(argc == 14) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
		checkParameter(argv[8], argv[9]);
		checkParameter(argv[10], argv[11]);
		checkParameter(argv[12], argv[13]);
	}
	if(argc == 16) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
		checkParameter(argv[8], argv[9]);
		checkParameter(argv[10], argv[11]);
		checkParameter(argv[12], argv[13]);
		checkParameter(argv[14], argv[15]);
	}
	if(argc == 18) {
		checkParameter(argv[4], argv[5]);
		checkParameter(argv[6], argv[7]);
		checkParameter(argv[8], argv[9]);
		checkParameter(argv[10], argv[11]);
		checkParameter(argv[12], argv[13]);
		checkParameter(argv[14], argv[15]);
		checkParameter(argv[16], argv[17]);
	}
	return 1;
}

bool createDirectories() {
	bool read_index = true;
	struct stat sb;
	if (stat(index_location.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		FILE *pFile;
		FILE *pFile2;
		FILE *pFile3;
		string build = SSTR(build_number);
		string sizes_location = index_location + "//sizes" + build;
		string sites_location = index_location + "//sites" + build;
		string info_location = index_location + "//index" + build + ".info";
		pFile = fopen( sizes_location.c_str() , "rb" );
		pFile2 = fopen(sites_location.c_str(), "rb");
		pFile3 = fopen(info_location.c_str(), "rb");

		if(!pFile || !pFile2 || !pFile3) {
			cout << "Index does not exist at location '" << index_location << "'." << endl;
			cout << "Index will be created at location '" << index_location << "'." << endl;
			read_index = false;
		}
		if(pFile) {
			fclose(pFile);
		}
		if(pFile2) {
			fclose(pFile2);
		}
		if(pFile3) {
			fclose(pFile3);
		}
	}
	else {
		cout << "Directory '" << index_location << "' does not exist." << endl;
		read_index = false;
		const string out = "mkdir -p " + index_location;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << index_location << "'." << endl;
			exit(-1);
		}
		else {
			cout << "Directory '" << index_location << "' created." << endl;
			cout << "Index will be stored in '" << index_location << "'." << endl;
		}
	}

	if(overwrite_index) {
		read_index = false;
	}
	return read_index;

}

void secondAlign(int **res,  map<int, FastaRead> &read_names, map<int, Result> &correct_results,
		string &whole_genome, vector<Read> &reads, bool read_index, bool add, vector<Read> &unaligned) {
	int *sizes = res[1];
	int *sites = res[3];

	length_of_sizes = res[0][0];
	length_of_sites = res[2][0];

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	cout << "Number of reads: "<< reads.size() << "\n";
	string build = SSTR(build_number);
	//unsigned int cutoff = whole_genome.size();

	pthread_t threads[thread_num];
	int rc;
	pthread_attr_t attr;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	vector<vector<Result> > set_of_results;
	for(int i = 0; i < thread_num; i++) {
		vector<Result> tmp;
		set_of_results.push_back(tmp);
	}
	vector<vector<Read> > unaligned_reads;
	for(int i = 0; i < thread_num; i++) {
		vector<Read> tmp;
		unaligned_reads.push_back(tmp);
	}
	ThreadData3 datas3[thread_num];
	int difference = reads.size()/thread_num;

	for(int i = 0; i < thread_num; i++) {
		int start =  i*difference;
		int stop =  i*difference + difference;
		datas3[i] = ThreadData3(i, sizes, sites, &reads, start, stop, &whole_genome, &set_of_results[i], &unaligned_reads[i]);
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

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Reads aligned: " << timefinal << " seconds." << endl;

	vector<Result> tmp_results;
	for(unsigned int i = 0; i < set_of_results.size(); i++) {
		for(unsigned int j = 0; j < set_of_results[i].size(); j++) {
			tmp_results.push_back(set_of_results[i][j]);
		}
	}
	set_of_results.clear();

	for(unsigned int i = 0; i < unaligned_reads.size(); i++) {
		for(unsigned int j = 0; j < unaligned_reads[i].size(); j++) {
			unaligned.push_back(unaligned_reads[i][j]);
		}
	}
	set_of_results.clear();

	cout << "Aligned reads: " << tmp_results.size() << endl;

	string addon = "x";
	if(add) {
		addon = "";
	}

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	writeResults(tmp_results, correct_results, outdir + "//" + "results" + addon + build + ".txt", !add, whole_genome);
	writeSamResults(tmp_results, reads, read_names, outdir + "//" + "results" + addon + build + ".sam");
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Results writen: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	Statistic s = Statistic(reads.size());
	writeStatistics(tmp_results, correct_results, outdir + "//" + "statistics" + addon + build + ".txt", s);
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

	/*if(first) {
		gettimeofday(&t1, NULL);
		startday = t1.tv_sec;
		startday2 = t1.tv_usec;
		vector<int> coverage;
		for(unsigned int i = 0; i < cutoff; i++) {
			coverage.push_back(0);
		}
		calculateCoverageFromResults(coverage, outdir + "//" + "results" + build + ".txt");
		writeCoverage(outdir + "//" + "coverage" + build + ".txt", coverage);
		writeCoverageFiltered(outdir + "//" + "coverage" + build + ".cut", coverage);
		gettimeofday(&t2, NULL);
		endday = t2.tv_sec;
		endday2 = t2.tv_usec;
		timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		gettimeofday(&t2, NULL);
		cout << "Coverage written: " << timefinal << " seconds. " << endl;

		gettimeofday(&t1, NULL);
		startday = t1.tv_sec;
		startday2 = t1.tv_usec;
		string part_genome;
		string in1 = outdir + "//" + "newref" + build + ".info";
		string in2 = outdir + "//" + "newref" + build + ".fa";
		writeNewRef(in1, in2, coverage, part_genome, whole_genome);
		whole_genome.clear();
		gettimeofday(&t2, NULL);
		endday = t2.tv_sec;
		endday2 = t2.tv_usec;
		timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "New ref written: " << timefinal << " seconds. " << endl;

		cout << "Creating second index..." << endl;
		gettimeofday(&t1, NULL);
		startday = t1.tv_sec;
		startday2 = t1.tv_usec;
		cout << "New ref written: " << timefinal << " seconds. " << endl;
		res = createIndex(true, part_genome, true, in2, index_location, 6, pow(2, 2*6), build_number);
		gettimeofday(&t2, NULL);
		endday = t2.tv_sec;
		endday2 = t2.tv_usec;
		timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;

		cout << "Size of new ref: " << part_genome.size() << endl;

		KEYLEN = 6;

		vector<Read> tmp_unaligned_reads;

		readUnalignedReads(tmp_unaligned_reads, reads, outdir + "//" + "results" + build + ".txt");

		cout << "Unaligned reads: " << tmp_unaligned_reads.size() << endl;

		secondAlign(res, read_names, correct_results, whole_genome, tmp_unaligned_reads);

		cout << "Program ended." << endl;
		free(sizes);
		free(sites);

		delete [] res[0];
		delete [] res[2];
		delete [] res;
	}*/
}

int main(int argc, char *argv[]) {

	if(readParams(argc, argv) < 0) exit(-1);

	if(precise) {
		MIN_QSCORE_MULT2 = 0.1;
		PRESCAN_QSCORE_THRESH = 0.57;
		Z_SCORE_MULT = 20;
		MIN_QSCORE_MULT = 0.025;
		MIN_SCORE_MULT = 0.15;
		DYNAMIC_QSCORE_THRESH_PERFECT = 0.8;
		DYNAMIC_QSCORE_THRESH = 0.6;
		DYNAMIC_SCORE_THRESH = 0.84;
	}

	string genome_ref = argv[3];
	string whole_genome;
	bool read_index = createDirectories();
	vector<Read> tmp_unaligned_reads;

	int **res;
	if(read_index) {
		cout << "Reading index..." << endl;
		res = readIndex(whole_genome, genome_ref, index_location, false);
	}
	else {
		cout << "Creating index..." << endl;
		res = createIndex(true, whole_genome, false, genome_ref, index_location, KEYLEN, keyspace, build_number);
	}

	cout << "genome size: " << whole_genome.size() << endl;

	vector<Read> reads;
	map<int, FastaRead> read_names;
	map<int, Result> correct_results;
	readReads(reads, read_names, correct_results, argv[2]);

	secondAlign(res, read_names, correct_results, whole_genome, reads, read_index, true, tmp_unaligned_reads);

	string build = SSTR(build_number);
	unsigned int cutoff = whole_genome.size();

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	vector<int> coverage;
	for(unsigned int i = 0; i < cutoff; i++) {
		coverage.push_back(0);
	}
	calculateCoverageFromResults(coverage, outdir + "//" + "results" + build + ".txt");
	writeCoverage(outdir + "//" + "coverage" + build + ".txt", coverage);
	writeCoverageFiltered(outdir + "//" + "coverage" + build + ".cut", coverage);
	gettimeofday(&t2, NULL);
	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;
	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	gettimeofday(&t2, NULL);
	cout << "Coverage written: " << timefinal << " seconds. " << endl;

	drugi = true;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	string part_genome;
	string in1 = outdir + "//" + "newref" + build + ".info";
	string in2 = outdir + "//" + "newref" + build + ".fa";
	writeNewRef(in1, in2, coverage, part_genome, whole_genome, positions_index);
	whole_genome.clear();
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "New ref written: " << timefinal << " seconds. " << endl;

	cout << "part_genome: " << part_genome.size() << endl;

	KEYLEN = 13;
	cout << "Creating second index..." << endl;
	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	res = createIndex(true, part_genome, true, in2, index_location, KEYLEN, pow(2, 2*KEYLEN), build_number);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;

	cout << "Size of new ref: " << part_genome.size() << endl;

	//readUnalignedReads(tmp_unaligned_reads, reads, outdir + "//" + "results" + build + ".txt");

	cout << "Unaligned reads: " << tmp_unaligned_reads.size() << endl;

	total_base_num = 0;
	for(unsigned int i = 0; i < tmp_unaligned_reads.size(); i++) {
		total_base_num += tmp_unaligned_reads[i].content.size();
	}

	aligned_base_num = 0;

	vector<Read> un_re;

/*
	MIN_QSCORE_MULT2 = 0.005;
	PRESCAN_QSCORE_THRESH = 0.6 * 0.95;
	Z_SCORE_MULT = 25;
	MIN_QSCORE_MULT = 0.005;
	MIN_SCORE_MULT = 0.02;
	DYNAMIC_SCORE_THRESH = 0.64;
	*/

	secondAlign(res, read_names, correct_results, part_genome, tmp_unaligned_reads, false, false, un_re);

	cout << positions_index[11590581] << endl;
	cout << positions_index[11618588] << endl;

	cout << "pos index: " << positions_index.size() << endl;

	string ovo =  outdir + "//" + "posint" + build + ".txt";
	ofstream posin(ovo.c_str());
	for(unsigned int i = 0; i < positions_index.size(); i++) {
		posin << positions_index[i] << " ";
		if(i % 100 == 0) posin << endl;
	}

	cout << "Program ended." << endl;

	return 0;
}

