# bfMEM

bfMEM is a tool for Maximal Exact Matches (MEMs) detecting. It is based on Bloom filter and rolling hash. The program is written in C++11 (tested with g++ >= 6.2.1) and works on Linux.


## Download & Compile

	git clone https://github.com/yuansliu/bfmem.git
	cd bfmem
	make

## Usage

	./bfmem -r H.all.fa -q M.all.fa -o hm-100.txt [options]

Parameters

	-r  	reference genome, a multi-FASTA file
	-q  	reference genome, a multi-FASTA file
	-o  	output file
	-l  	minimal length of matches; default is 100
	-k  	length of k-mer
	-t  	number of threads
	-s  	strands; 
			default is foward; 
			'r' is reverse-complement; 
			'b' is both foward and reverse-complement.
	-h  	print help message

## Status
Under review

## Citation
Yuansheng Liu, Leo Yu Zhang, Jinyan Li; bfMEM: fast detection of maximal exact matches via Bloom filter and rolling hash. 2018.

## Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>
