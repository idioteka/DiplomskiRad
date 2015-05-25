
#include "headers.h"
#include "config.h"

bool overlap(int a1, int b1, int a2, int b2){
	return a2<=b1 && b2>=a1;
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

void calculateStatistics3(Result &r1, Result &r2, Statistic &statistics) {

	int startOf = r1.start-r2.start;
	int stopOf = r1.stop-r2.stop;

	if(abs(r1.start-r2.start) < 15) {
		statistics.start15++;
	}
	if(abs(r1.stop-r2.stop) < 15) {
		statistics.stop15++;
	}
	if(abs(r1.start-r2.start) < 15 && abs(r1.stop-r2.stop) < 15) {
		statistics.start_and_stop15++;
	}
	if(abs(r1.start-r2.start) < 15 || abs(r1.stop-r2.stop) < 15) {
		statistics.start_or_stop15++;
	}
	if(abs(r1.start-r2.start) < 30) {
		statistics.start30++;
	}
	if(abs(r1.stop-r2.stop) < 30) {
		statistics.stop30++;
	}
	if(abs(r1.start-r2.start) < 30 && abs(r1.stop-r2.stop) < 30) {
		statistics.start_and_stop30++;
	}
	if(abs(r1.start-r2.start) < 30 || abs(r1.stop-r2.stop) < 30) {
		statistics.start_or_stop30++;
	}
	if(abs(r1.start-r2.start) < 50) {
		statistics.start50++;
	}
	if(abs(r1.stop-r2.stop) < 50) {
		statistics.stop50++;
	}
	if(abs(r1.start-r2.start) < 50 && abs(r1.stop-r2.stop) < 50) {
		statistics.start_and_stop50++;
	}
	if(abs(r1.start-r2.start) < 50 || abs(r1.stop-r2.stop) < 50) {
		statistics.start_or_stop50++;
	}
	if(abs(r1.start-r2.start) < 100) {
		statistics.start100++;
	}
	if(abs(r1.stop-r2.stop) < 100) {
		statistics.stop100++;
	}
	if(abs(r1.start-r2.start) < 100 && abs(r1.stop-r2.stop) < 100) {
		statistics.start_and_stop100++;
	}
	if(abs(r1.start-r2.start) < 100 || abs(r1.stop-r2.stop) < 100) {
		statistics.start_or_stop100++;
	}

}

void calculateStatistics5(Config &config, Result &r1, Result &r2, Statistic &statistics) {

	int startOf = config.POSITIONS_INDEX[r1.start]-r2.start;
	int stopOf = config.POSITIONS_INDEX[r1.stop]-r2.stop;

	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 15) {
		statistics.start15++;
	}
	if(abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 15) {
		statistics.stop15++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 15 && abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 15) {
		statistics.start_and_stop15++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 15 || abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 15) {
		statistics.start_or_stop15++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 30) {
		statistics.start30++;
	}
	if(abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 30) {
		statistics.stop30++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 30 && abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 30) {
		statistics.start_and_stop30++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 30 || abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 30) {
		statistics.start_or_stop30++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 50) {
		statistics.start50++;
	}
	if(abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 50) {
		statistics.stop50++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 50 && abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 50) {
		statistics.start_and_stop50++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 50 || abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 50) {
		statistics.start_or_stop50++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 100) {
		statistics.start100++;
	}
	if(abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 100) {
		statistics.stop100++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 100 && abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 100) {
		statistics.start_and_stop100++;
	}
	if(abs(config.POSITIONS_INDEX[r1.start]-r2.start) < 100 || abs(config.POSITIONS_INDEX[r1.stop]-r2.stop) < 100) {
		statistics.start_or_stop100++;
	}
}

void calculateStatistics4(Config &config, Result &r1, Result &r2, Statistic &statistics) {

	statistics.aligned_reads++;
	int sum = 0;

	if(r1.gapArray.size() == 0) {
		r1.gapArray.push_back(r1.start);
		r1.gapArray.push_back(r1.stop);
	}
	statistics.found_exons += r1.gapArray.size() / 2;
	statistics.total_exons += r2.gapArray.size() / 2;
	for(unsigned int i = 0; i < r1.gapArray.size(); i+=2) {
		bool overlaped = false;
		for(unsigned int j = 0; j < r2.gapArray.size(); j+=2) {
			if(overlap(r2.gapArray[j], r2.gapArray[j+1], config.POSITIONS_INDEX[r1.gapArray[i]], config.POSITIONS_INDEX[r1.gapArray[i+1]])) {
				int tmp = calculateOverlap(r2.gapArray[j], r2.gapArray[j+1], config.POSITIONS_INDEX[r1.gapArray[i]], config.POSITIONS_INDEX[r1.gapArray[i+1]]);
				sum += tmp;
				overlaped = true;
			}
		}
		if(overlaped) {
			statistics.overlaping_exons++;
		}
		statistics.finded_bases += (config.POSITIONS_INDEX[r1.gapArray[i+1]] - config.POSITIONS_INDEX[r1.gapArray[i]]);
	}
	statistics.covered_bases += sum;
}


void calculateStatistics2(Result &r1, Result &r2, Statistic &statistics) {

	statistics.aligned_reads++;

	int sum = 0;

	if(r1.gapArray.size() == 0) {
		r1.gapArray.push_back(r1.start);
		r1.gapArray.push_back(r1.stop);
	}
	statistics.found_exons += r1.gapArray.size() / 2;
	statistics.total_exons += r2.gapArray.size() / 2;
	for(unsigned int i = 0; i < r1.gapArray.size(); i+=2) {
		int num_over = 0;
		bool overlaped = false;
		for(unsigned int j = 0; j < r2.gapArray.size(); j+=2) {
			if(overlap(r2.gapArray[j], r2.gapArray[j+1], r1.gapArray[i], r1.gapArray[i+1])) {
				int tmp = calculateOverlap(r2.gapArray[j], r2.gapArray[j+1], r1.gapArray[i], r1.gapArray[i+1]);
				sum += tmp;
				num_over++;
				overlaped = true;

			}
		}
		if(overlaped) {
			statistics.overlaping_exons++;
		}
		statistics.finded_bases += (r1.gapArray[i+1] - r1.gapArray[i]);
		//if(num_over > 1) cout << "num over: " << num_over << endl;
	}
	statistics.covered_bases += sum;

}

void writeTotalStatistics(Config &config, string infile, string sssfile, Statistic &statistic) {
	ofstream out_stat(infile.c_str());

	out_stat << "Total number of reads: " << statistic.total_reads << endl;
	out_stat << "Number of aligned reads: " << statistic.aligned_reads << " - ";
	out_stat << (100 * statistic.aligned_reads / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Unaligned number of reads: " << statistic.unaligned_reads << " - ";
	out_stat << (100 * statistic.unaligned_reads / (double) statistic.total_reads) << "%" << endl;

	out_stat <<	"Total number of exons: " << statistic.total_exons << endl;
	out_stat << "Number of found exons: " << statistic.found_exons << " - ";
	out_stat << (100 * (statistic.found_exons / (double) statistic.total_exons)) << "%" << endl;
	out_stat << "Number of overlaping exons: " << statistic.overlaping_exons << " - ";
	out_stat << (100 * (statistic.overlaping_exons / (double) statistic.found_exons)) << "%" << endl;

	out_stat << "Reads with correct starts (<15): " << statistic.start15 << " - ";
	out_stat << (100 * statistic.start15 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start15 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct stops (<15): " << statistic.stop15 << " - ";
	out_stat << (100 * statistic.stop15 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.stop15 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts and stops (<15): " << statistic.start_and_stop15 << " - ";
	out_stat << (100 * statistic.start_and_stop15 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start_and_stop15 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts or stops (<15): " << statistic.start_or_stop15 << " - ";
	if(config.SPLIT_COUNT != 1) {
		out_stat << (100 * statistic.start_or_stop15 / (double) (statistic.aligned_reads * (config.SPLIT_COUNT/2))) << "% - ";
		out_stat << (100 * statistic.start_or_stop15 / (double) (statistic.total_reads * (config.SPLIT_COUNT/2))) << "%" << endl;
	}
	else {
		out_stat << (100 * statistic.start_or_stop15 / (double) (statistic.aligned_reads)) << "% - ";
		out_stat << (100 * statistic.start_or_stop15 / (double) (statistic.total_reads)) << "%" << endl;
	}
	out_stat << "Reads with correct starts (<30): " << statistic.start30 << " - ";
	out_stat << (100 * statistic.start30 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start30 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct stops (<30): " << statistic.stop30 << " - ";
	out_stat << (100 * statistic.stop30 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.stop30 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts and stops (<30): " << statistic.start_and_stop30 << " - ";
	out_stat << (100 * statistic.start_and_stop30 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start_and_stop30 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts or stops (<30): " << statistic.start_or_stop30 << " - ";
	if(config.SPLIT_COUNT != 1) {
		out_stat << (100 * statistic.start_or_stop30 / (double) (statistic.aligned_reads * (config.SPLIT_COUNT/2))) << "% - ";
		out_stat << (100 * statistic.start_or_stop30 / (double) (statistic.total_reads * (config.SPLIT_COUNT/2))) << "%" << endl;
	}
	else {
		out_stat << (100 * statistic.start_or_stop30 / (double) (statistic.aligned_reads)) << "% - ";
		out_stat << (100 * statistic.start_or_stop30 / (double) (statistic.total_reads)) << "%" << endl;
	}
	out_stat << "Reads with correct starts (<50): " << statistic.start50 << " - ";
	out_stat << (100 * statistic.start50 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start50 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct stops (<50): " << statistic.stop50 << " - ";
	out_stat << (100 * statistic.stop50 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.stop50 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts and stops (<50): " << statistic.start_and_stop50 << " - ";
	out_stat << (100 * statistic.start_and_stop50 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start_and_stop50 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts or stops (<50): " << statistic.start_or_stop50 << " - ";
	if(config.SPLIT_COUNT != 1) {
		out_stat << (100 * statistic.start_or_stop50 / (double) (statistic.aligned_reads * (config.SPLIT_COUNT/2))) << "% - ";
		out_stat << (100 * statistic.start_or_stop50 / (double) (statistic.total_reads * (config.SPLIT_COUNT/2))) << "%" << endl;
	}
	else {
		out_stat << (100 * statistic.start_or_stop50 / (double) (statistic.aligned_reads)) << "% - ";
		out_stat << (100 * statistic.start_or_stop50 / (double) (statistic.total_reads)) << "%" << endl;
	}
	out_stat << "Reads with correct starts (<100): " << statistic.start100 << " - ";
	out_stat << (100 * statistic.start100 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start100 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct stops (<100): " << statistic.stop100 << " - ";
	out_stat << (100 * statistic.stop100 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.stop100 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts and stops (<100): " << statistic.start_and_stop100 << " - ";
	out_stat << (100 * statistic.start_and_stop100 / (double) statistic.aligned_reads) << "% - ";
	out_stat << (100 * statistic.start_and_stop100 / (double) statistic.total_reads) << "%" << endl;
	out_stat << "Reads with correct starts or stops (<100): " << statistic.start_or_stop100 << " - ";
	if(config.SPLIT_COUNT != 1) {
		out_stat << (100 * statistic.start_or_stop100 / (double) (statistic.aligned_reads * (config.SPLIT_COUNT/2))) << "% - ";
		out_stat << (100 * statistic.start_or_stop100 / (double) (statistic.total_reads * (config.SPLIT_COUNT/2))) << "%" << endl;
	}
	else {
		out_stat << (100 * statistic.start_or_stop100 / (double) (statistic.aligned_reads)) << "% - ";
		out_stat << (100 * statistic.start_or_stop100 / (double) (statistic.total_reads)) << "%" << endl;
	}
	out_stat << "Total number of bases: " << statistic.total_bases << endl;
	out_stat << "Number of covered bases: " << statistic.covered_bases << endl;
	out_stat << "Number of founded bases: " << statistic.finded_bases << endl;
	out_stat << "Percentage of covered bases of all reads: ";
	out_stat << (100 * (statistic.covered_bases / (double) statistic.total_bases)) << "%" << endl;
	out_stat << "Percentage of covered bases of aligned reads: ";
	out_stat << (100 * (statistic.covered_bases / (double) statistic.aligned_bases)) << "%" << endl;

	out_stat << "Percentage of covered bases of founded reads: ";
	out_stat << (100 * (statistic.covered_bases / (double) statistic.finded_bases)) << "%" << endl;

	out_stat.close();
}

void addTotalStatistics(Config &config, vector<Result> &results, map<int, Result> &correct_results, Statistic &statistic) {

	for(unsigned int i = 0; i < results.size(); i++) {
		Result r1 = results[i];
		map<int, Result>::iterator pos = correct_results.find(r1.br);
		Result r2 = pos->second;

		if(config.IS_SECOND_PHASE) {
			calculateStatistics4(config, r1, r2, statistic);
			//cout << "4" << endl;
			calculateStatistics5(config, r1, r2, statistic);
		//	cout << "5" << endl;
		} else {
			calculateStatistics2(r1, r2, statistic);
			calculateStatistics3(r1, r2, statistic);
		}
	}
	if(config.IS_SECOND_PHASE) {
		statistic.aligned_reads = statistic.aligned_reads / (double) config.SPLIT_COUNT;
	}
	else {
		statistic.aligned_reads = statistic.aligned_reads / (double) config.SPLIT_COUNT;
	}

	statistic.unaligned_reads = statistic.total_reads - statistic.aligned_reads;
	//cout << "sss: " << statistic.unaligned_reads << endl;

}

void writeStatistics(Config &config, vector<Result> &results, map<int, Result> &correct_results, string infile, string sssfile, Statistic &statistic) {
	//ofstream out_stat(infile.c_str());

	statistic.total_bases = config.TOTAL_BASE_NUM;
	statistic.aligned_bases = config.ALIGNED_BASE_NUM;

	map<int, int> sss;
	map<int, int> sts;

	for(unsigned int i = 0; i < results.size(); i++) {
		Result r1 = results[i];
		map<int, Result>::iterator pos = correct_results.find(r1.br);
		Result r2 = pos->second;

		if(config.IS_SECOND_PHASE) {
			calculateStatistics4(config, r1, r2, statistic);
			calculateStatistics5(config, r1, r2, statistic);
		} else {
			calculateStatistics2(r1, r2, statistic);
			calculateStatistics3(r1, r2, statistic);
		}
	}
	if(config.IS_SECOND_PHASE) {
		statistic.aligned_reads = statistic.aligned_reads / (double) config.SPLIT_COUNT;
	}
	else {
		statistic.aligned_reads = statistic.aligned_reads / (double) config.SPLIT_COUNT;
	}
	statistic.unaligned_reads = statistic.total_reads - statistic.aligned_reads;

	writeTotalStatistics(config, infile, sssfile, statistic);

}

