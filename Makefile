run:
	g++ -g -I /home/donmeza2/software/htslib-1.20_install/include -D__STDC_LIMIT_MACROS -L /home/donmeza2/software/htslib-1.20_install/lib -std=c++11 -lhts -lz -lgcov -lpthread main.cpp Node.cpp Graph.cpp util.cpp -o main 
	./main /data/KolmogorovLab/donmeza2/cancer_phasing_data/chr10/HG002_HG00733_chr10_merged_140-160.vcf.gz /data/KolmogorovLab/donmeza2/cancer_phasing_data/chr10/HG002_HG00733_chr10_merged_140-160.bam
clean:
	rm main