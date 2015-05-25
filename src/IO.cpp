
#include "headers.h"
#include "config.h"
#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void writeReads(vector<Read> &reads, string dat) {
	ofstream out(dat.c_str());
	for (int var = 0; var < reads.size(); ++var) {
		Read r = reads[var];
		out << r.br << "\t" << r.content << "\n";
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void readUnalignedReads(vector <Read> &tmp_unaligned_reads, vector <Read> &reads, string infile) {
	ifstream ifs(infile.c_str());

	ofstream out("test.txt");

	string line;
	string line2;
	string line3;
	string line4;
	vector<int> brs;

	while(getline(ifs, line)) {
		getline(ifs, line2);
		getline(ifs, line3);
		getline(ifs, line4);

		char delim2 = '-';
		stringstream ss(line);
		string item;
		int br;
		while (std::getline(ss, item, delim2)) {
			br = atoi(item.c_str());
			out << br << endl;;
			break;
		}
		brs.push_back(br);
	}
	cout << "done" << endl;

	for(unsigned int i = 1; i <= reads.size(); i++) {
		if(std::find(brs.begin(), brs.end(), i) != brs.end()) {
			/* v contains x */
		} else {
			Read tm = reads[i-1];
			tmp_unaligned_reads.push_back(tm);
		}
	}
}

void gapArrayFromString(vector<long> &gapArray, string str) {
	char delim = ' ';
	char delim2 = '-';
	stringstream ss(str);
	string item;
	while (std::getline(ss, item, delim)) {
		if(item[item.size()-1] == ',') {
			item.erase(item.size()-1);
		}
		stringstream ss2(item);
		string first;
		string second;
		getline(ss2, first, delim2);
		getline(ss2, second, delim2);
		gapArray.push_back(atoi(first.c_str())-1);
		gapArray.push_back(atoi(second.c_str())-1);
	}
}

void readFastaReads(vector<Read> &reads, string infile, map<int, FastaRead> &fastaReads) {

	FILE *inReads = fopen(infile.c_str(), "r");
	if (inReads == NULL) {
		cout << "Fasta file does not exist: " << infile << endl;
		exit(-1);
	}
	char buff[256];
	string read;
	int br = 1;
	fgets(buff, 255, inReads);
	string idn = buff;
	idn.resize(idn.size() - 1);

	while (fgets(buff, 255, inReads)) {
		// ignore lines that start with genome desc, they start with '>'
		if (buff[0] != '>') {
			string tmp = buff;
			tmp.resize(tmp.size() - 1);  // remove endl
			read += tmp;
		} else {
			Read r(br, read);
			reads.push_back(r);
			string tmp = buff;
			tmp.resize(tmp.size() - 1);  // remove endl
			string first = idn;
			string second = read;
			FastaRead fr(first, second);
			idn = tmp;
			std::map<int,FastaRead>::iterator it = fastaReads.begin();
			fastaReads.insert(it, pair<int,FastaRead>(br, fr));
			br++;
			read.clear();
		}
	}

	Read r(br, read);
	reads.push_back(r);
	string first = idn;
	string second = read;
	FastaRead fr(first, second);
	std::map<int,FastaRead>::iterator it = fastaReads.begin();
	fastaReads.insert(it, pair<int,FastaRead>(br, fr));

	int test2 = fastaReads.size();
}


void readFastaReads2(vector<Read> &reads, string infile, map<int, FastaRead> fastaReads) {
	ifstream ifs(infile.c_str());
	if(!ifs) {
		cout << "File " << infile << " does not exist." << endl;
		exit(-1);
	}
	string line;
	string line2;
	int br = 1;
	while(getline(ifs, line)) {
		getline(ifs, line2);
		Read r(br, line2);
		reads.push_back(r);
		string first = line.substr(2, line.size()-2);
		string second = line2;
		FastaRead fr(first, second);
		std::map<int,FastaRead>::iterator it = fastaReads.begin();
		fastaReads.insert(it, pair<int,FastaRead>(br, fr));
		br++;
	}

}

void readReads(Config &config, vector<Read> &reads, map<int, FastaRead> &fastaReads, map<int, Result> &results, string infile) {
	ifstream ifs(infile.c_str());
	if(!ifs) {
		cout << "File " << infile << " does not exist." << endl;
		exit(-1);
	}
	string line;

	int br = 1;
	while(getline(ifs, line)) {
		vector<string> lines = split(line, '\t');
		string first = SSTR(br);
		string second = lines[5];
		FastaRead fr(first, second);
		std::map<int,FastaRead>::iterator it2 = fastaReads.begin();
		fastaReads.insert(it2, pair<int,FastaRead>(br, fr));
		Result r = Result(br, atoi(lines[1].c_str()), atoi(lines[2].c_str()), 0, 0);
		r.matchString = lines[3];
		gapArrayFromString(r.gapArray, lines[4]);
		std::map<int,Result>::iterator it = results.begin();
		results.insert(it, pair<int,Result>(br, r));
		Read rea(br, lines[5]);
		reads.push_back(rea);
		config.TOTAL_BASE_NUM += lines[5].size();
		br++;
	}
}


string compressMatchString(string &matchString) {
	char current = matchString[0];
	int number = 1;
	string compressed;
	for(unsigned int i = 1; i < matchString.size(); i++) {
		if(current == matchString[i]) {
			number++;
		}
		else if((current == 'm' && matchString[i] == 's') || (current == 's' && matchString[i] == 'm')) {
			number++;
		}
		else {
			if(current == 's') current = 'm';
			compressed += SSTR(number);
			compressed.push_back(current);
			number = 1;
			current = matchString[i];
		}
	}
	if(current == 's') current = 'm';
	compressed += SSTR(number);
	compressed.push_back(current);
	return compressed;
}

void sortResults(vector<Result> &results) {
	for(unsigned int i = 0; i < results.size()-1; i++) {
		unsigned int min = i;
		for(unsigned int j = i+1; j < results.size(); j++) {
			if(results[j].br < results[min].br) min = j;
		}
		if(min != i) {
			Result r = results[i];
			results[i] = results[min];
			results[min] = r;
		}
	}
}

string getSegmentName(int start, vector<ReferenceSegment> &referenceSegments) {
	string name;
	for(int i = referenceSegments.size()-1; i >= 0; i--) {
		ReferenceSegment rs = referenceSegments[i];
		if(rs.start < start) {
			return rs.name;
		}
	}
	return name;
}

void writeSamResults(vector<Result> &set_of_results, vector<Read> &reads, map<int, FastaRead> &read_names, string infile, vector<ReferenceSegment> &referenceSegments) {
	vector<Result> results;
	for(unsigned int j = 0; j < set_of_results.size(); j++) {
		results.push_back(set_of_results[j]);
	}
	sortResults(results);
	ofstream out_res(infile.c_str());
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		if(r.gapArray.size() == 0) {
			r.gapArray.push_back(r.start);
			r.gapArray.push_back(r.stop);
		}
		map<int, FastaRead>::iterator pos = read_names.find(r.br);
		FastaRead fr = pos->second;

		out_res << (fr.name) << "\t";
		for(unsigned j = 0; j < r.gapArray.size(); j+=2) {
			out_res << r.gapArray[j] << "-" << r.gapArray[j+1];
			if(j+2 != r.gapArray.size()) out_res << ",";
		}
		out_res << "\t";
		string segmentName = getSegmentName(r.start, referenceSegments);
		out_res << segmentName << "\t";
		out_res << r.start << "\t";
		out_res << r.stop << "\t";
		out_res << compressMatchString(r.matchString) << "\t";
		out_res << r.score << "\t";
		out_res << "\n";
	}
}

void writeResults(Config &config, vector<Result> &set_of_results, map<int, Result> &correct_results, string infile) {

	vector<Result> results;
	for(unsigned int j = 0; j < set_of_results.size(); j++) {
		results.push_back(set_of_results[j]);
	}
	cout << "size of results: " << results.size() << endl;
	sortResults(results);
	ofstream out_res(infile.c_str());
	int br = 1;
	for(unsigned int i = 0; i < results.size(); i++) {
		Result r = results[i];
		map<int, Result>::iterator pos = correct_results.find(r.br);
		Result r2 = pos->second;
		int start = r.start;
		int stop = r.stop;
		if(config.IS_SECOND_PHASE) {
			start = config.POSITIONS_INDEX[start];
			stop = config.POSITIONS_INDEX[stop];
			out_res << r.br << "-" << start << "-" << stop << "-" << (stop-start) << endl;
		}
		else {
			out_res << r.br << "-" << r.start << "-" << r.stop << "-" << (r.stop - r.start) << endl;
		}
		out_res << r2.br << "-" << r2.start << "-" << r2.stop << "-" << (r2.stop-r2.start) << endl;
		if(r.gapArray.size() == 0) {
			r.gapArray.push_back(r.start);
			r.gapArray.push_back(r.stop);
		}
		if(config.IS_SECOND_PHASE) {
			for(unsigned int j = 0; j < r.gapArray.size(); j+=2) {
				out_res << config.POSITIONS_INDEX[r.gapArray[j]] << "-" << config.POSITIONS_INDEX[r.gapArray[j+1]] << ", ";
			}
		}
		else {
			for(unsigned int j = 0; j < r.gapArray.size(); j+=2) {
				out_res << r.gapArray[j] << "-" << r.gapArray[j+1] << ", ";
			}
		}
		out_res << endl;
		out_res << r.matchString << endl;
		//out_res << whole_genome.substr(r.start, (r.stop - r.start)) << endl;
		br++;
	}

	out_res.close();
}


void writeInfo(Info &info, Config &config) {
	string filename = config.OUTDIR + "//" + "align.info";
	ofstream out_info(filename.c_str());

	out_info << "Elapsed execution time: " << info.total_execution_time << " seconds." << endl;
	out_info << "Index generation time: " << info.index_generation_time << " seconds." << endl;
	out_info << "First phase time: " << info.first_align_time2 << " seconds." << endl;
	out_info << "Coverage time: " << info.coverage_time << " seconds." << endl;
	out_info << "New reference time: " << info.new_ref_time << " seconds." << endl;
	out_info << "Second ref indexing time: " << info.second_ref_indexing_time << " seconds." << endl;
	out_info << "Second align time: " << info.second_align_time2 << " seconds." << endl;
	out_info << "Third ref indexing time: " << info.third_ref_indexing_time << " seconds." << endl;
	out_info << "Third align time: " << info.third_align_time2 << " seconds." << endl;

	out_info << "Thread number: " << info.number_of_threads << endl;
	out_info << "Whole genome size: " << info.genome_size << " bases." << endl;
	out_info << "Total number of reads: " << info.read_number << endl;
	out_info << "Number of aligned reads: " << info.aligned_read_number << endl;
	out_info << "First phase data: " << endl;
	out_info << "\t - keylen: " << info.first_phase_keylen << endl;
	out_info << "\t - max indel: " << info.first_phase_max_indel << endl;
	out_info << "\t - precision: " << info.first_phase_multy_precision << endl;
	out_info << "\t - read length: " << info.first_phase_readlen << endl;
	out_info << "Second phase data: " << endl;
	out_info << "\t - keylen: " << info.second_phase_keylen << endl;
	out_info << "\t - max indel: " << info.second_phase_max_indel << endl;
	out_info << "\t - precision: " << info.second_phase_multy_precision << endl;
	out_info << "\t - read length: " << info.second_phase_readlen << endl;
	out_info << "Third phase data: " << endl;
	out_info << "\t - keylen: " << info.third_phase_keylen << endl;
	out_info << "\t - max indel: " << info.third_phase_max_indel << endl;
	out_info << "\t - precision: " << info.third_phase_multy_precision << endl;
	out_info << "\t - read length: " << info.first_phase_readlen << endl;
	out_info << "Mode: " << info.second_phase_mode << endl;
	out_info << "Coverage gaplen: " << info.coverage_gaplen << endl;
	out_info << "Coverage cutoff: " << info.coverage_threshold << endl;
	out_info << "Coverage padding: " << info.coverage_padding << endl;
	out_info << "New reference size: " << info.new_ref_size << endl;
	out_info << "Unprecised reads: " << endl;
	out_info << "\t - better: " << info.better << endl;
	out_info << "\t - worse: " << info.worse << endl;
	out_info << "\t - same: " << info.same << endl;
	out_info << "\t - unfound: " << info.unfound << endl;
	out_info << "Memory after first align: " << info.memory_after_first_align << endl;
	out_info << "Memory after second align: " << info.memory_after_second_align << endl;
	out_info << "Memory after third align: " << info.memory_after_third_align << endl;
	out_info << "Memory at the end: " << info.memory_at_the_end << endl;
}

