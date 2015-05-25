
#include "config.h"

bool overlap(int a1, int b1, int a2, int b2);
int calculateGappedPrecisePercentage(int precise_start, int precise_stop, int start1, int stop1, int start2, int stop2, string &matchString);
int calculateOverlap(int start1, int stop1, int start2, int stop2);
void calculateStatistics3(Result &r1, Result &r2, Statistic &statistics, map<int, int> &sss, map<int, int> &sts);
void calculateStatistics5(Config &config, Result &r1, Result &r2, Statistic &statistics, map<int, int> &sss, map<int, int> &sts);
void calculateStatistics4(Config &config, Result &r1, Result &r2, Statistic &statistics);
void calculateStatistics2(Result &r1, Result &r2, Statistic &statistics);
void writeTotalStatistics(Config &config, string infile, string sssfile, Statistic &statistic);
void addTotalStatistics(Config &config, vector<Result> &results, map<int, Result> &correct_results, Statistic &statistic);
void writeStatistics(Config &config, vector<Result> &results, map<int, Result> &correct_results, string infile, string sssfile, Statistic &statistic);
