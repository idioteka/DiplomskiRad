
#include "config.h"

bool overlap2(int a1, int b1, int a2, int b2);
int calculateOverlap2(int start1, int stop1, int start2, int stop2);
void calculateStatistics32(Result &r1, Result &r2, Statistic &statistics, map<int, int> &sss, map<int, int> &sts);
void calculateStatistics22(Result &r1, Result &r2, Statistic &statistics);
void writeTotalStatistics2(Config &config, string infile, string sssfile, Statistic &statistic);
void writeStatistics2(Config &config, vector<Result> &results, map<int, Result> &correct_results, string infile, string sssfile, Statistic &statistic);
