
#include "config.h"

int key_num;
static int KEYLEN;
static int KEYSPACE;
static int BUILD_NUMBER;

struct ThreadData {
	int thread_id;
	long *array;
	string *whole_genome;
	ThreadData() {
	}
	ThreadData(int thread_id_, long *array_, string *whole_genome_) {
		thread_id = thread_id_;
		array = array_;
		whole_genome = whole_genome_;
	}
};

struct ThreadData2 {
	int thread_id;
	long *array;
	long *sizes;
	string *whole_genome;
	ThreadData2() {

	}
	ThreadData2(int thread_id_, long *array_, long *sizes_, string &whole_genome_) {
		thread_id = thread_id_;
		array = array_;
		sizes = sizes_;
		whole_genome = &whole_genome_;
	}
};

void extractGenomeFromFile(string genome_file, string &whole_genome, vector<ReferenceSegment> &referenceSegments) {

	//whole_genome = "";
	FILE *inGenome = fopen(genome_file.c_str(), "r");
	if (inGenome == NULL) {
		cout << "Reference file does not exist: " << genome_file << endl;
		exit(-1);
	}
	char buff[256];

	int oldPosition = 0;
	int position = 0;
	string name;

	fgets(buff, 255, inGenome);
	string tmp = buff;
	tmp.resize(tmp.size() - 1);
	name = tmp;
	int br = 1;

	while (fgets(buff, 255, inGenome)) {
		// ignore lines that start with genome desc, they start with '>'
		if (buff[0] != '>') {
			string tmp = buff;
			tmp.resize(tmp.size() - 1);  // remove endl
			whole_genome += tmp;
			position += tmp.size();
		} else {
			ReferenceSegment rs = ReferenceSegment(br, name, oldPosition);
			oldPosition = position;
			referenceSegments.push_back(rs);
			br++;

			string tmp = buff;
			tmp.resize(tmp.size() - 1);
			name = tmp;
		}
	}

	ReferenceSegment rs = ReferenceSegment(br, name, oldPosition);
	referenceSegments.push_back(rs);
}

int getKeyFromKmer(string &whole_genome, long start, long stop) {
	int key = 0;
	for(long i = start; i < stop; i++) {
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

	long start = 0;
	long sum = 0;

	while(getCodeFromBase(td->whole_genome[0][start]) < 0) {
		start++;
	}

	for(long i = start; i < td->whole_genome[0].size()-(KEYLEN-1); i++) {
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
	long start = 0;

	while(getCodeFromBase(td->whole_genome[0][start]) < 0) {
		start++;
	}

	for(long i = start; i < td->whole_genome[0].size()-(KEYLEN-1); i++) {
		if(getCodeFromBase(td->whole_genome[0][i]) == (*td).thread_id) {
			int code = getCodeFromBase(td->whole_genome[0][i]);
			if(code >= 0) {
				int key = getKeyFromKmer(td->whole_genome[0], i, i+KEYLEN);
				if(key >= 0) {
					long location = (*td).sizes[key];
					(*td).array[location] = i;
					(*td).sizes[key]++;
				}
			}
		}
	}
	pthread_exit(NULL);
}

void writeSizes(long *sizes, string loc) {
	FILE* pFile;
	pFile = fopen(loc.c_str(), "wb");
	//TODO check if file exist
	fwrite(sizes, sizeof(long), KEYSPACE, pFile);
	fclose(pFile);
}

void writeSites(long *sites, long sum, string loc) {
	FILE* pFile;
	pFile = fopen(loc.c_str(), "wb");
	//TODO check if file exist
	fwrite(sites, sizeof(long), sum, pFile);
	fclose(pFile);
}

long *readArray(string filename, bool write_sum) {
	FILE * pFile;
	long lSize;
	long * buffer;

	pFile = fopen ( filename.c_str() , "rb" );

	if(!pFile) {
		cout << "File " << filename << " does not exit. Please create index." << endl;
		exit(-1);
	}

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	if(write_sum) {
		key_num = lSize/8;
	}
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (long*) malloc (sizeof(char)*lSize);

	// copy the file into the buffer:
	fread (buffer,1,lSize,pFile);
	/* the whole file is now loaded in the memory buffer. */

	// terminate
	fclose (pFile);
	return buffer;
}

long ** readIndex(string &whole_genome, string genome_ref, string index_loc, bool part_genome, Config &config, vector<ReferenceSegment> &referenceSegments) {

	KEYLEN = config.KEYLEN;
	KEYSPACE = config.KEYSPACE;
	BUILD_NUMBER = config.BUILD_NUMBER;

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;
	extractGenomeFromFile(genome_ref, whole_genome, referenceSegments);
	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;
	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Genome read: " << timefinal << " seconds." << endl;

	t1, t2;
	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;

	string build = SSTR(BUILD_NUMBER);
	string loc = index_loc + "//" + "sizes" + build;
	if(part_genome) {
		loc = index_loc + "//" + "6sizes" + build;
	}
	long *sizes = readArray(loc, false);

	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Sizes loaded: " << timefinal << " seconds." << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	loc = index_loc + "//" + "sites" + build;
	if(part_genome) {
		loc = index_loc + "//" + "6sites" + build;
	}
	long *sites = readArray(loc, true);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;

	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Sites loaded: " << timefinal << " seconds." << endl;

	long **result = new long*[4];
	result[0] = new long[1];
	result[0][0] = KEYSPACE;
	result[1] = sizes;
	result[2] = new long[1];
	result[2][0] = key_num;
	result[3] = sites;
	return result;
}

void writeInfo(long sizes, long sites, string path, string &genome) {
	//string build = SSTR(build_number);
	ofstream ofs(path.c_str());
	ofs << "build number: " << BUILD_NUMBER << endl;
	ofs << "genome size: " << genome.size() << endl;
	ofs << "sizes size: " << sizes << endl;
	ofs << "sites size: " << sites << endl;
	ofs << "KEYLEN: " << KEYLEN << endl;
}

long ** createIndex(bool write_to_file, string &whole_genome, bool part_genome, string genome_ref, string index_loc, Config &config, vector<ReferenceSegment> &referenceSegments) {

	KEYLEN = config.KEYLEN;
	KEYSPACE = config.KEYSPACE;
	BUILD_NUMBER = config.BUILD_NUMBER;

	if(!part_genome) {
		extractGenomeFromFile(genome_ref, whole_genome, referenceSegments);
	}

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	pthread_t threads[4];
	int rc;
	long *sizes = new long[KEYSPACE] {0};
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

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;

	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Created sizes: " << timefinal << endl;


	long sum=0;
	for(long i = 0; i < KEYSPACE; i++) {
		long temp = sizes[i];
		sizes[i] = sum;
		sum += temp;
	}

	if(write_to_file) {
		timeval t1, t2;
		gettimeofday(&t1, NULL);
		long startday = t1.tv_sec;
		long startday2 = t1.tv_usec;
		string build = SSTR(BUILD_NUMBER);
		string loc = index_loc + "//sizes" + build;
		if(part_genome) {
			loc = index_loc + "//6sizes" + build;
		}
		writeSizes(sizes, loc);
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

	cout << "sum: " << sum << endl;

	long *sites = new long[sum];
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

	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;
	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "Created sites: " << timefinal << endl;

	string build = SSTR(BUILD_NUMBER);
	if(write_to_file) {
		gettimeofday(&t1, NULL);
		long startday = t1.tv_sec;
		long startday2 = t1.tv_usec;
		string loc = index_loc + "//sites" + build;
		if(part_genome) {
			loc = index_loc + "//6sites" + build;
		}
		writeSites(sites, sum, loc);
		gettimeofday(&t2, NULL);
		long endday = t2.tv_sec;
		long endday2 = t2.tv_usec;
		double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "Wrote sites: " << timefinal << endl;
	}

	for(long i = KEYSPACE-1; i > 0; i--) {
		sizes[i] = sizes[i-1];
	}
	sizes[0] = 0;

	string loc = index_loc + "//index" + build + ".info";
	if(part_genome) {
		loc = index_loc + "//6index" + build + ".info";
	}
	writeInfo(KEYSPACE, sum, loc, whole_genome);

	long **result = new long*[4];
	result[0] = new long[1];
	result[0][0] = KEYSPACE;
	result[1] = sizes;
	result[2] = new long[1];
	result[2][0] = sum;
	result[3] = sites;
	return result;

}

int main4(int argc, char *argv[]) {
	return 0;
}


