#include "headers.h"
#include "config.h"

/*
 *
 * OFFSETS
 *
 */

int getDesiredKeyNumber(Config &config, int readlen, float density) {
	int slots = readlen - config.KEYLEN +1;
	int desired = (int) ceil(( readlen * density) / config.KEYLEN);
	desired = max(config.MIN_KEYS_DESIRED, desired);
	desired = min(slots, desired);
	return desired;
}

void makeOffsets(Config &config, string read, vector<int> &offsets) {
	int test = read.size();
	float keyDen2 = (( config.MAX_DESIRED_KEYS * config.KEYLEN ) / (float) read.size());
	keyDen2 = max(config.MIN_KEY_DENSITY, keyDen2);
	keyDen2 = min(config.KEY_DENSITY, keyDen2);
	float keyDen3;
	if(read.size() <= 50){
		keyDen3 = config.MAX_KEY_DENSITY;
	}else if(read.size() >= 200){
		keyDen3 = config.MAX_KEY_DENSITY-0.5f;
	}else{
		keyDen3 = config.MAX_KEY_DENSITY - 0.003333333333f * (read.size()-50);
	}
	keyDen3 = max(config.KEY_DENSITY, keyDen3);

	int desiredKeysNumber = getDesiredKeyNumber(config, read.size(), keyDen2);

	float interval = (read.size() - config.KEYLEN) / (float) (max(desiredKeysNumber - 1, 1));
	float f = 0;

	for(int i=0, j = 0; i < desiredKeysNumber; i++){
		if(((unsigned int) round(f+interval)) > read.size()) {
			break;
		}
		offsets.push_back(j);

		f += interval;
		j = min(((int) read.size() - config.KEYLEN), (max(j+1, (int) round(f))));
	}
}

/*
 *
 * KEY CREATION
 *
 */


int getKeyFromKmer(Config &config, long start, long stop, string &read) {
	int key = 0;
	for(long i = start; i < stop; i++) {
		int code = getCodeFromBase(read[i]);
		if(code < 0) {
			return -1;
		}
		key = ((key<<2) | code);
	}
	return key;
}

void getReadKeys(Config &config, string read, vector<int> &offsets, vector<int> &keys) {
	for(unsigned int i = 0; i < offsets.size(); i++) {
		keys.push_back(getKeyFromKmer(config, offsets[i], offsets[i] + config.KEYLEN, read));
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

int countKeyHits(Config &config, int key, long *sizes) {
	int x;
	if(key+1 == config.LENGTH_OF_SIZES) {
		x = config.LENGTH_OF_SITES - sizes[key];
	}
	else {
		x = sizes[key+1] - sizes[key];
	}
	int rkey = reverseComplementBinary(key, config.KEYLEN);
	if(key == rkey) return x;
	else {
		int y;
		if(rkey+1 == config.LENGTH_OF_SIZES) {
			y = config.LENGTH_OF_SITES - sizes[rkey];
		}
		else {
			y = sizes[rkey+1] - sizes[rkey];
		}
		return x + y;
	}
}

int countHits(Config &config, vector<int> &keys, long *sizes, int max_len) {
	vector<int> keys_temp = keys;
	int num_hits = 0;
	for(unsigned int i = 0; i < keys_temp.size(); i++){
		int key = keys_temp[i];
		if(key >= 0){
			int len = countKeyHits(config, key, sizes);
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

int trimHits(Config &config, vector<int> &keys, vector<int> &keys_temp, long *sizes, int num_hits) {
	if(num_hits > 0) {
		int trigger = (3 * keys_temp.size())/4;
		if(num_hits < 20 && num_hits < trigger) {
			for(unsigned int i = 0; i < keys_temp.size(); i++) {
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(config, keys_temp, sizes, (config.MAX_LEN * 3) / 2);
		}
		if(num_hits < 18 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(config, keys_temp, sizes, config.MAX_LEN*2);
		}
		if(num_hits < 16 && num_hits < trigger){
			for(unsigned int i=0; i< keys.size(); i++){
				keys_temp[i] = keys[i];
			}
			num_hits = countHits(config, keys_temp, sizes, config.MAX_LEN*3);
		}
		if(num_hits < 14 && num_hits < trigger){
			for(unsigned int i=0; i<keys.size(); i++){
				keys_temp[i] = keys[i];}
			num_hits = countHits(config, keys_temp, sizes, config.MAX_LEN*5);
		}
	}
	return num_hits;
}

int getHits(Config &config, vector<int> &keys, int max_len, vector<long> &starts, vector<long> &stops, long *sizes){
	int num_hits = 0;
	for(unsigned int i = 0; i < keys.size(); i++) {
		int key = keys[i];
		starts.push_back(-1);
		stops.push_back(-1);
		if(key >= 0){
			int len = countKeyHits(config, key, sizes);
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

