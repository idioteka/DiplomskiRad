
int ** createIndex(bool write_to_file, string &whole_genome, string genome_ref, string index_loc, int keylen, int kspace, int build_num);
int ** readIndex(string &whole_genome, string genome_ref, string index_loc);
int getCodeFromBase(char base);
