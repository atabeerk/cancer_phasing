run:
	g++ -std=c++17 -g \
	-I/home/donmeza2/software/htslib-1.20_install/include \
	-Iinclude \
	-D__STDC_LIMIT_MACROS \
	-L/home/donmeza2/software/htslib-1.20_install/lib \
	main.cpp Node.cpp Graph.cpp connectSNVs.cpp util.cpp config.cpp \
	-lhts -lz -lpthread -lstdc++fs \
	-o main

clean:
	rm main