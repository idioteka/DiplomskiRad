
#include "headers.h"
#include "config.h"

/*
 *
 * PARAMS FUNCTIONS
 *
 */

void checkParameter(Info &info, Config &config, string command, string value) {
	if(strcmp(command.c_str(), "-special") == 0) {
		config.PRECISION_CUTOFF = atof(value.c_str());
		cout << "Precision cutoff set to: " << config.PRECISION_CUTOFF << endl;
	}
	else if(strcmp(command.c_str(), "-t") == 0) {
		config.THREAD_NUM = atoi(value.c_str());
		cout << "Number of threads set to " << config.THREAD_NUM << "." << endl;
		info.number_of_threads = config.THREAD_NUM;
	}
	else if(strcmp(command.c_str(), "-i") == 0) {
		config.INDEX_LOCATION = value;
		cout << "Index location set to '" << config.INDEX_LOCATION << "'." << endl;
	}
	else if(strcmp(command.c_str(), "-b") == 0) {
		config.BUILD_NUMBER = atoi(value.c_str());
		cout << "Build number set to " << value << "." << endl;
	}
	else if(strcmp(command.c_str(), "-c") == 0) {
		config.OVERWRITE_INDEX = true;
		config.BUILD_NUMBER = atoi(value.c_str());
		cout << "Index will be created. If index was previously created at ";
		cout << "set location it will be overwriten." << endl;
		cout << "Build number overwriten to " << value << "." << endl;
	}
	else if(strcmp(command.c_str(), "-s") == 0) {
		config.SPLIT_COUNT = atoi(value.c_str());
		cout << "Reads will be split at " << value << " parts." << endl;
	}
	else if(strcmp(command.c_str(), "-ss") == 0) {
		config.SECOND_PHASE_SPLIT_COUNT = atoi(value.c_str());
		cout << "In second phase reads will further be split at " << value << " parts." << endl;
	}
	else if(strcmp(command.c_str(), "-k") == 0) {
		config.KEYLEN = atoi(value.c_str());
		config.KEYSPACE = pow(2, 2*config.KEYLEN);
		cout << "KEYLEN set to " << value << "." << endl;
		info.first_phase_keylen = config.KEYLEN;
	}
	else if(strcmp(command.c_str(), "-k2") == 0) {
		config.KEYLEN2 = atoi(value.c_str());
		cout << "KEYLEN2 set to " << value << "." << endl;
		info.second_phase_keylen = config.KEYLEN2;
	}
	else if(strcmp(command.c_str(), "-k3") == 0) {
		config.KEYLEN3 = atoi(value.c_str());
		cout << "KEYLEN3 set to " << value << "." << endl;
		info.third_phase_keylen = config.KEYLEN3;
	}
	else if(strcmp(command.c_str(), "-p") == 0) {
		config.PRECISE = false;
		config.MULTY_PRECISION = atoi(value.c_str());
		cout << "Set to low precision." << endl;
		info.first_phase_precision = config.PRECISE;
		info.first_phase_multy_precision = config.MULTY_PRECISION;
	}
	else if(strcmp(command.c_str(), "-sp") == 0) {
		config.SECOND_PHASE_PRECISE = false;
		config.SECOND_PHASE_MULTY_PRECISION = atoi(value.c_str());
		cout << "Second phase set to low precision." << endl;
		info.second_phase_precision = config.SECOND_PHASE_PRECISE;
		info.second_phase_multy_precision = config.SECOND_PHASE_MULTY_PRECISION;
	}
	else if(strcmp(command.c_str(), "-tp") == 0) {
		config.THIRD_PHASE_MULTY_PRECISION = atoi(value.c_str());
		cout << "Third phase set to low precision." << endl;
		info.third_phase_multy_precision = config.THIRD_PHASE_MULTY_PRECISION;
	}
	else if(strcmp(command.c_str(), "-mi") == 0) {
		config.MAX_INDEL2 = atoi(value.c_str());
		cout << "Max indel set to " << value << "." << endl;
		info.first_phase_max_indel = config.MAX_INDEL2;
	}
	else if(strcmp(command.c_str(), "-smi") == 0) {
		config.SECOND_PHASE_MAX_INDEL2 = atoi(value.c_str());
		cout << "Second phase max indel set to " << value << "." << endl;
		info.second_phase_max_indel = config.SECOND_PHASE_MAX_INDEL2;
	}
	else if(strcmp(command.c_str(), "-tmi") == 0) {
		config.THIRD_PHASE_MAX_INDEL2 = atoi(value.c_str());
		cout << "Third phase max indel set to " << value << "." << endl;
		info.third_phase_max_indel = config.THIRD_PHASE_MAX_INDEL2;
	}
	else if(strcmp(command.c_str(), "-ct") == 0) {
		config.COV_THRES = atoi(value.c_str());
		cout << "Coverage threshold set to " << value << "." << endl;
		info.coverage_threshold = config.COV_THRES;
	}
	else if(strcmp(command.c_str(), "-cp") == 0) {
		config.COV_PADDING = atoi(value.c_str());
		cout << "Coverage padding set to " << value << "." << endl;
		info.coverage_padding = config.COV_PADDING;
	}
	else if(strcmp(command.c_str(), "-cg") == 0) {
		config.COV_GAPLEN = atoi(value.c_str());
		cout << "Coverage gaplen set to " << value << "." << endl;
		info.coverage_gaplen = config.COV_GAPLEN;
	}
	else if(strcmp(command.c_str(), "-spm") == 0) {
		config.SECOND_PHASE_MODE = atoi(value.c_str());
		cout << "Second phase mode set to " << value << "." << endl;
		info.second_phase_mode = config.SECOND_PHASE_MODE;
	}
	else if(strcmp(command.c_str(), "-pc") == 0) {
		config.PRECISION_CUTOFF = atof(value.c_str());
		cout << "Precision cutoff set to " << value << "." << endl;
	}
	else {
		cout << "Parameter " << command << " not recognized." << endl;
		exit(-1);
	}
}

int readParams(int argc, char *argv[], Config &config, Info &info) {
	if(argc < 4) {
		cout << "Program must run with arguments: <destination_folder> <reads_file> <genome_reference_file> -t <thread_number> -i <index_location> -c <create_index>" << endl;
		cout << "<destination_folder> - folder where results will be stored." << endl;
		cout << "<reads_file> - file with reads in fasta format." << endl;
		cout << "<genome_reference_file> - file with regerence genome." << endl;
		cout << "-t <thread_number> - Optional parameter. Number of threads, default 4." << endl;
		cout << "-k <KEYLEN> - Optional parameter. Length of the key. Default 13." << endl;
		cout << "-b <build_number> - Optional parameter. Build number of the index. By default set to 1." << endl;
		cout << "-i <index_location> - Optional parameter. Location of created index if index exists at <index_location>. If index does not exist, program creates index at <index_location>. " << endl;
		cout << "By default program looks for index in the <destination_folder> location. If there is no index, by default program creates new index in <destination_folder>" << endl;
		cout << "If you wish to create index even if it exists use -c <create_index> option." << endl;
		cout << "-c <build_number> - Optional parameter. Use if you wish to create index even if index exist at given location. This will overwrite previously created index at given location. It should be called with build number; for example: -c 2." << endl;
		cout << "-s <split_count> - Optional parameter. At how many party every read should be split. Default no split." << endl;
		cout << "-p <precision> - Optimal parameter. By default set to more precise, if set reads are aligned less precise." << endl;
		cout << "" << endl;
		return -1;
	}
	struct stat sb;
	cout << "Program started: " << endl;

	if (stat(argv[1], &sb) == 0 && S_ISDIR(sb.st_mode)) {
		config.OUTDIR = argv[1];
		cout << "Results will be stored in '" << config.OUTDIR << "'." << endl;
	}
	else {
		config.OUTDIR = argv[1];
		const string out = "mkdir -p " + config.OUTDIR;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << config.OUTDIR << "'." << endl;
			return -1;
		}
		else {
			cout << "Directory '" << config.OUTDIR << "' created." << endl;
			cout << "Results will be stored in '" << config.OUTDIR << "'." << endl;
		}
	}

	config.INDEX_LOCATION = config.OUTDIR;

	for(int i = 0; i < (argc - 4)/2 ; i++) {
		checkParameter(info, config, argv[4 + i*2], argv[5 + i*2]);
	}
	return 1;
}

bool createDirectories(Config &config) {
	bool read_index = true;
	struct stat sb;
	if (stat(config.INDEX_LOCATION.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		FILE *pFile;
		FILE *pFile2;
		FILE *pFile3;
		string build = SSTR(config.BUILD_NUMBER);
		string sizes_location = config.INDEX_LOCATION + "//sizes" + build;
		string sites_location = config.INDEX_LOCATION + "//sites" + build;
		string info_location = config.INDEX_LOCATION + "//index" + build + ".info";
		pFile = fopen( sizes_location.c_str() , "rb" );
		pFile2 = fopen(sites_location.c_str(), "rb");
		pFile3 = fopen(info_location.c_str(), "rb");

		if(!pFile || !pFile2 || !pFile3) {
			cout << "Index does not exist at location '" << config.INDEX_LOCATION << "'." << endl;
			cout << "Index will be created at location '" << config.INDEX_LOCATION << "'." << endl;
			read_index = false;
		}
		if(pFile) {
			fclose(pFile);
		}
		if(pFile2) {
			fclose(pFile2);
		}
		if(pFile3) {
			fclose(pFile3);
		}
	}
	else {
		cout << "Directory '" << config.INDEX_LOCATION << "' does not exist." << endl;
		read_index = false;
		const string out = "mkdir -p " + config.INDEX_LOCATION;
		if(system(out.c_str())) {
			cout << "Unable to create directory '" << config.INDEX_LOCATION << "'." << endl;
			exit(-1);
		}
		else {
			cout << "Directory '" << config.INDEX_LOCATION << "' created." << endl;
			cout << "Index will be stored in '" << config.INDEX_LOCATION << "'." << endl;
		}
	}

	if(config.OVERWRITE_INDEX) {
		read_index = false;
	}
	return read_index;

}


void updatePrecision2(Config &config, int prec) {
	if(prec == 1) {
		config.MIN_QSCORE_MULT2 = 0.1;
		config.PRESCAN_QSCORE_THRESH = 0.57;
		config.Z_SCORE_MULT = 20;
		config.MIN_QSCORE_MULT = 0.025;
		config.MIN_SCORE_MULT = 0.15;
		config.DYNAMIC_SCORE_THRESH = 0.84;
		config.MAX_DESIRED_KEYS = 15;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 1.5;
		config.KEY_DENSITY = 1.9;
		config.MAX_KEY_DENSITY = 3.0;
		config.SLOW_ALIGN_PADDING=6;
		config.MAXIMUM_MAX_HITS_REDUCTION=5;
		config.POINTS_MATCH = 70;
		config.POINTS_MATCH2 = 100;
		config.POINTS_SUB = -127;
		config.POINTS_SUB2 = -51;
		config.POINTS_SUB3 = -25;
		config.POINTS_INS = -395;
		//MAX_INDEL = 100;
		//MAX_INDEL2 = 8*100;
	}
	else if(prec == 2) {
		config.MIN_QSCORE_MULT2 = 0.1;
		config.PRESCAN_QSCORE_THRESH = 0.57;
		config.Z_SCORE_MULT = 20;
		config.MIN_QSCORE_MULT = 0.025;
		config.MIN_SCORE_MULT = 0.15;
		config.DYNAMIC_SCORE_THRESH = 0.84;
		config.MAX_DESIRED_KEYS = 63;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 2.8;
		config.KEY_DENSITY = 3.5;
		config.MAX_KEY_DENSITY = 4.5;
		config.SLOW_ALIGN_PADDING=8;
		config.MAXIMUM_MAX_HITS_REDUCTION=6;
	}
	else if(prec == 3) {
		config.MIN_QSCORE_MULT2 = 0.005;
		config.PRESCAN_QSCORE_THRESH = 0.6 * 0.95;
		config.Z_SCORE_MULT = 25;
		config.MIN_QSCORE_MULT = 0.005;
		config.MIN_SCORE_MULT = 0.02;
		config.DYNAMIC_SCORE_THRESH = 0.64;
		config.MAX_DESIRED_KEYS = 15;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 1.5;
		config.KEY_DENSITY = 1.9;
		config.MAX_KEY_DENSITY = 3.0;
		config.SLOW_ALIGN_PADDING=6;
		config.MAXIMUM_MAX_HITS_REDUCTION=5;
	}
	else {
		config.MIN_QSCORE_MULT2 = 0.005;
		config.PRESCAN_QSCORE_THRESH = 0.6 * 0.95;
		config.Z_SCORE_MULT = 25;
		config.MIN_QSCORE_MULT = 0.005;
		config.MIN_SCORE_MULT = 0.02;
		config.DYNAMIC_SCORE_THRESH = 0.64;
		config.MAX_DESIRED_KEYS = 63;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 2.8;
		config.KEY_DENSITY = 3.5;
		config.MAX_KEY_DENSITY = 4.5;
		config.SLOW_ALIGN_PADDING=8;
		config.MAXIMUM_MAX_HITS_REDUCTION=6;
		config.POINTS_MATCH = 90;
		config.POINTS_MATCH2 = 100;
		config.POINTS_SUB = -137;
		config.POINTS_SUB2 = -49;
		config.POINTS_SUB3 = -25;
		config.POINTS_INS = -205;
		//MAX_INDEL = 100;
		//MAX_INDEL2 = 8*100;
	}
}

void updatePrecision(Config &config, bool prec) {
	if(prec) {
		config.MIN_QSCORE_MULT2 = 0.1;
		config.PRESCAN_QSCORE_THRESH = 0.57;
		config.Z_SCORE_MULT = 20;
		config.MIN_QSCORE_MULT = 0.025;
		config.MIN_SCORE_MULT = 0.15;
		config.DYNAMIC_SCORE_THRESH = 0.84;
		config.MAX_DESIRED_KEYS = 15;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 1.5;
		config.KEY_DENSITY = 1.9;
		config.MAX_KEY_DENSITY = 3.0;
		config.SLOW_ALIGN_PADDING=6;
		config.MAXIMUM_MAX_HITS_REDUCTION=5;
		config.POINTS_MATCH = 70;
		config.POINTS_MATCH2 = 100;
		config.POINTS_SUB = -127;
		config.POINTS_SUB2 = -51;
		config.POINTS_SUB3 = -25;
		config.POINTS_INS = -395;
		//MAX_INDEL = 100;
		//MAX_INDEL2 = 8*100;
	}
	else {
		config.MIN_QSCORE_MULT2 = 0.005;
		config.PRESCAN_QSCORE_THRESH = 0.6 * 0.95;
		config.Z_SCORE_MULT = 25;
		config.MIN_QSCORE_MULT = 0.005;
		config.MIN_SCORE_MULT = 0.02;
		config.DYNAMIC_SCORE_THRESH = 0.64;
		config.MAX_DESIRED_KEYS = 63;
		config.MIN_KEYS_DESIRED = 2;
		config.MIN_KEY_DENSITY = 2.8;
		config.KEY_DENSITY = 3.5;
		config.MAX_KEY_DENSITY = 4.5;
		config.SLOW_ALIGN_PADDING=8;
		config.MAXIMUM_MAX_HITS_REDUCTION=6;
		config.POINTS_MATCH = 90;
		config.POINTS_MATCH2 = 100;
		config.POINTS_SUB = -137;
		config.POINTS_SUB2 = -49;
		config.POINTS_SUB3 = -25;
		config.POINTS_INS = -205;
		//MAX_INDEL = 100;
		//MAX_INDEL2 = 8*100;
	}
}
