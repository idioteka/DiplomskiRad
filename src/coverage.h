

/*
 *
 * FILTER REF
 *
 */
void writeNewRef(Config &config, string file, string file2, vector<int> &coverage, string &part_genome, string &whole_genome, vector<int> &positions);

/*
 *
 * COVERAGE
 *
 */
void calculateCoverageFromCig(vector<int> &coverage, string infile);
void calculateCoverageFromResults(vector<int> &coverage, string infile);
void calculateCoverageFromVector(Config &config, vector<int> &coverage, vector<Result> &results);
void writeCoverageFiltered(string path, vector <int> &coverage);
void writeCoverage(string path, vector <int> &coverage);
void checkCoverage(vector<int> &coverage);

void evaluateCoverage(vector<int> &final_coverage, vector<int> &correct_coverage, string infile);
