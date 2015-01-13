
#include "config.h"

int key_num;

struct ThreadData {
	int thread_id;
	int *array;
	string *whole_genome;
	ThreadData() {
	}
	ThreadData(int thread_id_, int *array_, string *whole_genome_) {
		thread_id = thread_id_;
		array = array_;
		whole_genome = whole_genome_;
	}
};

struct ThreadData2 {
	int thread_id;
	int *array;
	int *sizes;
	string *whole_genome;
	ThreadData2() {

	}
	ThreadData2(int thread_id_, int *array_, int *sizes_, string &whole_genome_) {
		thread_id = thread_id_;
		array = array_;
		sizes = sizes_;
		whole_genome = &whole_genome_;
	}
};

void extractGenomeFromFile(string genome_file, string &whole_genome) {

	//whole_genome = "";
	FILE *inGenome = fopen(genome_file.c_str(), "r");
	if (inGenome == NULL) {
		cout << "Reference file does not exist: " << genome_file << endl;
		exit(-1);
	}
	char buff[256];

	while (fgets(buff, 255, inGenome)) {
		// ignore lines that start with genome desc, they start with '>'
		if (buff[0] != '>') {
			string tmp = buff;
			tmp.resize(tmp.size() - 1);  // remove endl
			whole_genome += tmp;
		}
	}
	//return whole_genome;
}

int getCodeFromBase(char base) {
	char lower_base = tolower(base);
	if(lower_base == 'a') {
		return 0;
	} else if(lower_base == 'c') {
		return 1;
	}
	else if(lower_base == 'g') {
		return 2;
	}
	else if(lower_base == 't') {
		return 3;
	}
	else return -1;
}

int getKeyFromKmer(string &whole_genome, int start, int stop) {
	int key = 0;
	for(int i = start; i < stop; i++) {
		int code = getCodeFromBase(whole_genome[i]);
		if(code < 0) {
			return -1;
		}
		key = ((key<<2) | code);
	}
	return key;
}


void *countKeys(void *threadid) {
	ThreadData *td = (ThreadData *) threadid;

	int start = 0;
	int sum = 0;

	while(getCodeFromBase(td->whole_genome[0][start]) < 0) {
		start++;
	}

	for(unsigned int i = start; i < td->whole_genome[0].size()-(KEYLEN-1); i++) {
		if(getCodeFromBase(td->whole_genome[0][i]) == (*td).thread_id) {
			int code = getCodeFromBase(td->whole_genome[0][i]);
			if(code >= 0) {
				int key = getKeyFromKmer(td->whole_genome[0], i, i+KEYLEN);
				if(key >= 0) {
					(*td).array[key]++;
					sum++;
				}
			}
		}
	}

	pthread_exit(NULL);
}


void *fillArrays(void *threadid) {
	ThreadData2 *td = (ThreadData2 *) threadid;
	int start = 0;

	while(getCodeFromBase(td->whole_genome[0][start]) < 0) {
		start++;
	}

	for(unsigned int i = start; i < td->whole_genome[0].size()-(KEYLEN-1); i++) {
		if(getCodeFromBase(td->whole_genome[0][i]) == (*td).thread_id) {
			int code = getCodeFromBase(td->whole_genome[0][i]);
			if(code >= 0) {
				int key = getKeyFromKmer(td->whole_genome[0], i, i+KEYLEN);
				if(key >= 0) {
					int location = (*td).sizes[key];
					(*td).array[location] = i;
					(*td).sizes[key]++;
				}
			}
		}
	}
	pthread_exit(NULL);
}

void writeSizes(int *sizes, string loc) {
	FILE* pFile;
	pFile = fopen(loc.c_str(), "wb");
	//TODO check if file exist
	fwrite(sizes, sizeof(int), keyspace, pFile);
	fclose(pFile);
}

void writeSites(int *sites, int sum, string loc) {
	FILE* pFile;
	pFile = fopen(loc.c_str(), "wb");
	//TODO check if file exist
	fwrite(sites, sizeof(int), sum, pFile);
	fclose(pFile);
}

int *readArray(string filename, bool write_sum) {
	FILE * pFile;
	long lSize;
	int * buffer;

	pFile = fopen ( filename.c_str() , "rb" );

	if(!pFile) {
		cout << "File " << filename << " does not exit. Please create index." << endl;
		exit(-1);
	}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	if(write_sum) {
		key_num = lSize/4;
	}
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (int*) malloc (sizeof(char)*lSize);

	//cout << "Sizes size: " << lSize << endl;

	// copy the file into the buffer:
	fread (buffer,1,lSize,pFile);
	/* the whole file is now loaded in the memory buffer. */

	// terminate
	fclose (pFile);
	return buffer;
}

int ** readIndex(string &whole_genome, string genome_ref, string index_loc) {

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;
	extractGenomeFromFile(genome_ref, whole_genome);
	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;
	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Genome read: " << timefinal << " seconds." << endl;

	t1, t2;
	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;

	string build = SSTR(build_number);
	int *sizes = readArray(index_loc + "//" + "sizes" + build, false);

	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Sizes loaded: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	int *sites = readArray(index_loc + "//" + "sites" + build, true);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;

	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Sites loaded: " << timefinal << " seconds." << endl;

	int **result = new int*[4];
	result[0] = new int[1];
	result[0][0] = keyspace;
	result[1] = sizes;
	result[2] = new int[1];
	result[2][0] = key_num;
	result[3] = sites;
	return result;
}

void writeInfo(int sizes, int sites, string path) {
	string build = SSTR(build_number);
	ofstream ofs((path + "//index" + build + ".info").c_str());
	ofs << "build number:" << build_number << endl;
	ofs << "sizes size:" << sizes << endl;
	ofs << "sites size:" << sites << endl;
	ofs << "KEYLEN:" << KEYLEN << endl;
}

int ** createIndex(bool write_to_file, string &whole_genome, string genome_ref, string index_loc, int keylen, int kspace, int build_num) {

	KEYLEN = keylen;
	keyspace = kspace;
	build_number = build_num;

	extractGenomeFromFile(genome_ref, whole_genome);

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	pthread_t threads[4];
	int rc;
	int *sizes = new int[keyspace] {0};
	pthread_attr_t attr;
	void *status;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	ThreadData threadDatas[4];

	for(int i = 0; i < 4; i++) {
		threadDatas[i] = ThreadData(i, sizes, &whole_genome);
	}

	for(int i = 0; i < 4; i++) {
		//ThreadData td = ThreadData(i, sizes);
		rc = pthread_create(&threads[i], NULL, countKeys, (void *) &threadDatas[i]);
	}
	pthread_attr_destroy(&attr);

	for(int i=0; i < 4; i++ ){
		rc = pthread_join(threads[i], &status);
		if (rc){
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id :" << i ;
		cout << "  exiting with status :" << status << endl;
	}

	//int *sizes = countKeys(whole_genome);

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Created sizes: " << timefinal << endl;

	int sum=0;
	for(int i = 0; i < keyspace; i++) {
		int temp = sizes[i];
		sizes[i] = sum;
		sum += temp;
	}

	if(write_to_file) {
		timeval t1, t2;
		gettimeofday(&t1, NULL);
		long startday = t1.tv_sec;
		long startday2 = t1.tv_usec;
		string build = SSTR(build_number);
		writeSizes(sizes, index_loc + "//sizes" + build);
		gettimeofday(&t2, NULL);
		long endday = t2.tv_sec;
		long endday2 = t2.tv_usec;

		double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "Wrote sizes: " << timefinal << endl;
	}



	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;


	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t threads2[4];
	ThreadData2 threadDatas2[4];
	int *sites = new int[sum];
	for(int i = 0; i < 4; i++) {
		threadDatas2[i] =ThreadData2(i, sites, sizes, whole_genome);
	}
	for(int i = 0; i < 4; i++) {
		rc = pthread_create(&threads2[i], NULL, fillArrays, (void *)&threadDatas2[i]);
	}
	pthread_attr_destroy(&attr);
	for(int i=0; i < 4; i++ ){
		rc = pthread_join(threads2[i], &status);
		if (rc){
			cout << "Error:unable to join," << rc << endl;
			exit(-1);
		}
		cout << "Main: completed thread id :" << i ;
		cout << "  exiting with status :" << status << endl;
	}

	//int *sites = fillArrays(sizes, sum, whole_genome);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Created sites: " << timefinal << endl;

	if(write_to_file) {
		gettimeofday(&t1, NULL);
		long startday = t1.tv_sec;
		long startday2 = t1.tv_usec;
		string build = SSTR(build_number);
		writeSites(sites, sum, index_loc + "//sites" + build);
		gettimeofday(&t2, NULL);
		long endday = t2.tv_sec;
		long endday2 = t2.tv_usec;
		double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "Wrote sites: " << timefinal << endl;
	}

	for(int i = keyspace-1; i > 0; i--) {
		sizes[i] = sizes[i-1];
	}
	sizes[0] = 0;

	writeInfo(keyspace, sum, index_loc);

	int **result = new int*[4];
	result[0] = new int[1];
	result[0][0] = keyspace;
	result[1] = sizes;
	result[2] = new int[1];
	result[2][0] = sum;
	result[3] = sites;
	return result;

}

int main4(int argc, char *argv[]) {
	return 0;
}


