run:
	g++ -std=c++17 -g \
	-I/home/donmeza2/software/htslib-1.20_install/include \
	-Iinclude \
	-D__STDC_LIMIT_MACROS \
	-L/home/donmeza2/software/htslib-1.20_install/lib \
	main.cpp snv_analysis.cpp process_snv_output.cpp \
	-lhts -lz -lpthread -lstdc++fs \
	-o main

test:
	./main /data/KolmogorovLab/donmeza2/cancer_phasing_data/1437_T_chr7_118m-120m_gdrive/H1437.haplotagged.chr7_118000000_120000000.merged.vcf.gz \
	/data/KolmogorovLab/donmeza2/cancer_phasing_data/1437_T_chr7_118m-120m_gdrive/H1437.haplotagged.chr7_118000000_120000000.merged.bam \
	/data/KolmogorovLab/donmeza2/cancer_phasing_data/output/H1437_chr7_118-120_output

clean:
	rm -f main \
	rm -r /data/KolmogorovLab/donmeza2/cancer_phasing_data/output/H1437_chr7_118-120_output