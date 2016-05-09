# RNA Aligner

A program for alligning RNA sequencing reads to reference.

## Installation

To compile and create an executable run:
    g++ ./src/*.cpp -lpthread -o RNAAligner

Run with:
    RNAAligner <destination_folder> <reads_file> <genome_reference_file> -t <thread_number> -i <index_location> -c <create_index>

    <destination_folder> - folder where results will be stored.
    <reads_file> - file with reads in fasta format.
    <genome_reference_file> - file with regerence genome.
    -t <thread_number> - Optional parameter. Number of threads, default 4.
    -k <KEYLEN> - Optional parameter. Length of the key. Default 13.
    -b <build_number> - Optional parameter. Build number of the index. By default set to 1.
    -i <index_location> - Optional parameter. Location of created index if index exists at <index_location>. If index does not exist, program creates index at <index_location>. 
    By default program looks for index in the <destination_folder> location. If there is no index, by default program creates new index in <destination_folder>
    If you wish to create index even if it exists use -c <create_index> option.
    -c <build_number> - Optional parameter. Use if you wish to create index even if index exist at given location. This will overwrite previously created index at given location. It should be called with build number; for example: -c 2.
    -s <split_count> - Optional parameter. At how many party every read should be split. Default no split.
    -p <precision> - Optimal parameter. By default set to more precise, if set reads are aligned less precise.
