# bfMEM

bfMEM is a tool for Maximal Exact Matches (MEMs) detection. It is based on Bloom filter and rolling hash. The program is written in C++11 (tested with g++ >= 6.2.1) and works on Linux.


## Download & Compile

	git clone https://github.com/yuansliu/bfmem.git
	cd bfmem
	make

## Usage

	./bfmem -r H.all.fa -q M.all.fa -o hm-100.txt [options]

Parameters

	-r  	reference genome, a multi-FASTA file
	-q  	query genome, a multi-FASTA file
	-o  	output file

Options

	-l  	minimal length of matches; default is 100
	-k  	length of k-mer (some default values are set)
	-t  	number of threads
	-s  	strands; 
			default is foward; 
			'r' is reverse-complement; 
			'b' is both foward and reverse-complement.
	-h  	print help message

## Status
Submitted

## Citation
Yuansheng Liu, Leo Yu Zhang, Jinyan Li. Fast detection of maximal exact matches via fixed sampling of query k-mers and Bloom filtering of index k-mers. 2019.

## Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>
