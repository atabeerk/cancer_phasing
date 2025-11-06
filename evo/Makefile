run:
	g++ -std=c++17 -g \
	-I/Users/donmeza2/Documents/cancer_phasing/include/htslib \
	-Iinclude \
	-D__STDC_LIMIT_MACROS \
	-L/Users/donmeza2/Documents/cancer_phasing/include/htslib/lib \
	snv_cooccur.cpp snv_analysis.cpp process_snv_output.cpp \
	../include/htslib/lib/libhts.a -lz -lbz2 -llzma -lcurl -lpthread \
	-o snv_occur

test:
	./snv_occur /Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437_T_chr7_118m-120m_gdrive/H1437.haplotagged.chr7_118000000_120000000.merged.vcf.gz \
	/Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437_T_chr7_118m-120m_gdrive/H1437.haplotagged.chr7_118000000_120000000.merged.bam \
	H1437_chr7_118-120_out

bd:
	./snv_occur /Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437_T_chr1_100m-250m/H1437.haplotagged.chr1_100000000_250000000.merged.vcf.gz \
	/Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437_T_chr1_100m-250m/H1437.haplotagged.chr1_100000000_250000000.merged.bam \
	H1437_chr1_100-150_out_is_error

1437_chr2:
	./snv_occur /Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437T_chr2/H1437.haplotagged.chr2_1_242193528.hp1_only.vcf.gz \
	/Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data/1437T_chr2/H1437.haplotagged.chr2_1_242193528.merged.bam \
	H1437_chr2