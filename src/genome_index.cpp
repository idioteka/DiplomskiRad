
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
	clock_t sat = clock();
	ThreadData *td = (ThreadData *) threadid;

	//int *sizes = (*td).array;
	//cout << "threadid: " << td->thread_id << "\n";
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
/*
int *countKeys(string &whole_genome) {
	int *sizes = new int[keyspace] {0};
	int start = 0;
	int sum = 0;

	while(getCodeFromBase(whole_genome[start]) < 0) {
		start++;
	}

	for(unsigned int i = start; i < whole_genome.size()-(KEYLEN-1); i++) {

		int code = getCodeFromBase(whole_genome[i]);
		if(code >= 0) {
			int key = getKeyFromKmer(whole_genome, i, i+KEYLEN);
			if(key >= 0) {
				sizes[key]++;
				sum++;
			}
		}
	}

	return sizes;
}
*/


void *fillArrays(void *threadid) {
	ThreadData2 *td = (ThreadData2 *) threadid;
	int start = 0;
	//cout << "fill threadid: " << td->thread_id << "\n";
	//int *sites = (*td).array;
	//int *sizes = (*td).array;

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
/*
int *fillArrays(int *sizes, int sum, string &whole_genome) {
	int *sites = new int[sum];
	int start = 0;

	while(getCodeFromBase(whole_genome[start]) < 0) {
		start++;
	}

	for(unsigned int i = start; i < whole_genome.size()-(KEYLEN-1); i++) {

		int code = getCodeFromBase(whole_genome[i]);
		if(code >= 0) {
			int key = getKeyFromKmer(whole_genome, i, i+KEYLEN);
			if(key >= 0) {
				int location = sizes[key];
				sites[location] = i;
				sizes[key]++;
			}
		}
	}
	return sites;
}
*/

void writeSizes(int *sizes) {
	FILE* pFile;
	pFile = fopen("sizes", "wb");
	fwrite(sizes, sizeof(int), keyspace, pFile);
	fclose(pFile);
}

void writeSizes2(int *sizes) {
	ofstream sizes_output_stream(sizes_file_name.c_str());
	sizes_output_stream << keyspace;
	sizes_output_stream << "\n";
	for(int i = 0; i < keyspace; i++) {
		sizes_output_stream << sizes[i];
		sizes_output_stream << "\n";
	}
}

void writeSites(int *sites, int sum) {
	FILE* pFile;
	pFile = fopen("sites", "wb");
	fwrite(sites, sizeof(int), sum, pFile);
	fclose(pFile);
}

void writeSites2(int *sites, int sum) {
	ofstream sizes_output_stream(sites_file_name.c_str());
	sizes_output_stream << sum;
	sizes_output_stream << "\n";
	for(int i = 0; i < sum; i++) {
		sizes_output_stream << sites[i];
		sizes_output_stream << "\n";
	}
}

int *readArray(string filename, bool write_sum) {
	FILE * pFile;
	long lSize;
	int * buffer;

	pFile = fopen ( filename.c_str() , "rb" );

	// obtain file size:
	fseek (pFile , 0 , SEEK_END);
	lSize = ftell (pFile);
	if(write_sum) {
		key_num = lSize/4;
	}
	rewind (pFile);

	// allocate memory to contain the whole file:
	buffer = (int*) malloc (sizeof(char)*lSize);

	// copy the file into the buffer:
	fread (buffer,1,lSize,pFile);
	/* the whole file is now loaded in the memory buffer. */

	// terminate
	fclose (pFile);
	return buffer;
}

int *readArray2(string filename, bool write_sum) {
	ifstream ifs(filename.c_str());
	int size;
	ifs >> size;
	if(write_sum) key_num = size;
	int *array = new int[size];
	for(int i = 0; i < size; i++) {
		ifs >> array[i];
	}
	return array;
}

int ** readIndex(string &whole_genome) {

	extractGenomeFromFile(genome_file_name, whole_genome);

	timeval t1, t2;
	gettimeofday(&t1, NULL);
	long startday = t1.tv_sec;
	long startday2 = t1.tv_usec;

	int *sizes = readArray("sizes", false);

	gettimeofday(&t2, NULL);
	long endday = t2.tv_sec;
	long endday2 = t2.tv_usec;
	double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "READ SIZES: " << timefinal << endl;

	gettimeofday(&t1, NULL);
	startday = t1.tv_sec;
	startday2 = t1.tv_usec;
	int *sites = readArray("sites", true);
	gettimeofday(&t2, NULL);
	endday = t2.tv_sec;
	endday2 = t2.tv_usec;

	timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
	cout << "READ SITES: " << timefinal << endl;

	//length_of_sizes = keyspace;
	//length_of_sites = key_num;

	//length_of_sites = key_num;
	//cout << "length of sizes: " << length_of_sizes << endl;
	//cout << "length of sites: " << length_of_sites << endl;

	int **result = new int*[4];
	result[0] = new int[1];
	result[0][0] = keyspace;
	result[1] = sizes;
	result[2] = new int[1];
	result[2][0] = key_num;
	result[3] = sites;
	return result;
}

int ** createIndex(bool write_to_file, string &whole_genome) {

	extractGenomeFromFile(genome_file_name, whole_genome);

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
	cout << "CREATED SIZES: " << timefinal << endl;

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
		writeSizes(sizes);
		gettimeofday(&t2, NULL);
		long endday = t2.tv_sec;
		long endday2 = t2.tv_usec;

		double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "WROTE SIZES: " << timefinal << endl;
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
	cout << "CREATED SITES: " << timefinal << endl;

	if(write_to_file) {
		gettimeofday(&t1, NULL);
		long startday = t1.tv_sec;
		long startday2 = t1.tv_usec;
		writeSites(sites, sum);
		gettimeofday(&t2, NULL);
		long endday = t2.tv_sec;
		long endday2 = t2.tv_usec;
		double timefinal = ((endday - startday) * 1000000.0 + (endday2 - startday2))/ 1000000;
		cout << "WROTE SITES: " << timefinal << endl;
	}

	for(int i = keyspace-1; i > 0; i--) {
		sizes[i] = sizes[i-1];
	}
	sizes[0] = 0;

	int **result = new int*[4];
	result[0] = new int[1];
	result[0][0] = keyspace;
	result[1] = sizes;
	result[2] = new int[1];
	result[2][0] = sum;
	result[3] = sites;
	return result;

	return result;
}

int main2(int argc, char *argv[]) {

	//cout << gi->id;

	//ifstream ifs(genome_file_name.c_str());
//	string line;
	//int br = 1;

//	whole_genome = extractGenomeFromFile(genome_file_name);

//	int *sizes = countKeys();
/*	clock_t begin_read_sizes = clock();

	int *sizes = readArray(sizes_file_name, false);

	clock_t end_read_sizes = clock();
	cout << "ELAPSED READ SIZES: " << (double)(begin_read_sizes - end_read_sizes) / CLOCKS_PER_SEC << endl;
*/
/*
	int tops[5] = {0};
	int tops_keys[5] = {0};

	for(int i = 0; i < keyspace; i++) {
		if(tops[4] < sizes[i] && tops[3] < sizes[i] &&
				tops[2] < sizes[i] && tops[1] < sizes[i] && tops[0] < sizes[i]) {
			tops[4] = tops[3];
			tops[3] = tops[2];
			tops[2] = tops[1];
			tops[1] = tops[0];
			tops[0] = sizes[i];
			tops_keys[4] = tops_keys[3];
			tops_keys[3] = tops_keys[2];
			tops_keys[2] = tops_keys[1];
			tops_keys[1] = tops_keys[0];
			tops_keys[0] = i;
		}
		else if(tops[4] < sizes[i] && tops[3] < sizes[i] &&
				tops[2] < sizes[i] && tops[1] < sizes[i]) {
			tops[4] = tops[3];
			tops[3] = tops[2];
			tops[2] = tops[1];
			tops[1] = sizes[i];
			tops_keys[4] = tops_keys[3];
			tops_keys[3] = tops_keys[2];
			tops_keys[2] = tops_keys[1];
			tops_keys[1] = i;
		}
		else if (tops[4] < sizes[i] && tops[3] < sizes[i] &&
				tops[2] < sizes[i]) {
			tops[4] = tops[3];
			tops[3] = tops[2];
			tops[2] = sizes[i];
			tops_keys[4] = tops_keys[3];
			tops_keys[3] = tops_keys[2];
			tops_keys[2] = i;
		}
		else if(tops[4] < sizes[i] && tops[3] < sizes[i]) {
			tops[4] = tops[3];
			tops[3] = sizes[i];
			tops_keys[4] = tops_keys[3];
			tops_keys[3] = i;
		}
		else if(tops[4] < sizes[i]) {
			tops[4] = sizes[i];
			tops_keys[4] = i;
		}
	}

	for(int i = 0; i < 5; i++) {
		cout << tops[i] << " " << tops_keys[i] << "\n";
	}
*/
/*
	int sum=0;
	for(int i = 0; i < keyspace; i++) {
		int temp = sizes[i];
		sizes[i] = sum;
		sum += temp;
	}
*/
	//writeSizes(sizes);

	//cout << "sum: " << sum;

	//int *sites = fillArrays(sizes, sum);

	/*clock_t begin_read_sites = clock();

	int *sites = readArray(sites_file_name, true);

	clock_t end_read_sites = clock();
	cout << "ELAPSED READ SITES: " << (double)(begin_read_sites - end_read_sites) / CLOCKS_PER_SEC << endl;
*/

//	int **res = readIndex();

	return 0;
	//cout << res[0][0] << " " << res[2][0] << "\n";
	//cout << res[1][10] << " " << res[3][10] << "\n";

//	cout << sizes[0] << " " << sites[0];

//	writeSites(sites, sum);

//	for(int i = 0; i < 10; i++) {
//		cout << sites[i] << "\n";
//	}

	//while(getline(ifs, line)) {
		//if(line.substr(0,1) == ">") {
		//	cout << br;
			//cout << line;
			//cout << "\n";
		//}
//		br++;
	//}
}


