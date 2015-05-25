#include "headers.h"
#include "config.h"

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

vector<vector<long> > shrink(vector<long> &starts, vector<long> &stops, vector<int> &offsets, vector<int> &keys, int len){
	int numHits=0;
	for(int i=0; i<len; i++){
		if(starts[i]>=0){numHits++;}
	}
	vector<vector<long> > res;

	if(numHits == (int)offsets.size()){
		return res;
	}else{
		vector<long> offsets2;
		vector<long> starts2;
		vector<long> stops2;

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
bool checkIfAllBasesCovered(Config &config, vector<int> &offsets, string &read) {
	bool all_bases_covered = true;
	if(offsets[0] != 0) {
		all_bases_covered = false;
	}
	else if(offsets[(int)offsets.size() - 1] != ((int)read.size() - config.KEYLEN)) {
		all_bases_covered = false;
	}
	else{
		for(unsigned int i = 1; i < offsets.size(); i++){
			if(offsets[i] > offsets[i-1] + config.KEYLEN){
				all_bases_covered = false;
				break;
			}
		}
	}
	return all_bases_covered;
}

bool calculatePretendAllBasesCovered(Config &config, vector<int> &read_keys,vector<int> &read_keys_final,
				vector<int> &offsets, bool all_bases_covered, string &read) {
	bool pretend_all_bases_covered =
		(all_bases_covered ||
		read_keys_final.size() >= read_keys.size()-4 ||
		(read_keys_final.size() >= 9 &&
				(offsets[offsets.size()-1]-offsets[0]+config.KEYLEN)>max(40, (int)(read.size()*.75f))));
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
 * PROCESSES
 *
 */
int parseLine(char* line){
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

int getValue(){ //Note: this value is in KB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = parseLine(line);
			break;
		}
	}
	fclose(file);
	return result;
}

