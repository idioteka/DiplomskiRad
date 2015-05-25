

void readUnalignedReads(vector <Read> &tmp_unaligned_reads, vector <Read> &reads, string infile);
void gapArrayFromString(vector<long> &gapArray, string str);
void readFastaReads(vector<Read> &reads, string infile, map<int, FastaRead> &fastaReads);
void readFastaReads2(vector<Read> &reads, string infile, map<int, FastaRead> fastaReads);
void readReads(Config &config, vector<Read> &reads, map<int, FastaRead> &fastaReads, map<int, Result> &results, string infile);

void sortResults(vector<Result> &results);
string compressMatchString(string &matchString);
void writeSamResults(vector<Result> &set_of_results, vector<Read> &reads, map<int, FastaRead> &read_names, string infile, vector<ReferenceSegment> &referenceSegments);
void writeResults(Config &config, vector<Result> &set_of_results, map<int, Result> &correct_results, string infile);

void writeInfo(Info &info, Config &config);
void writeReads(vector<Read> &reads, string dat);
