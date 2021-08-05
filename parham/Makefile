default : libparham.so

lib%.so : %.cc
	g++ -g -O3 -fopenmp $< -shared -fPIC -D_GLIBCXX_PARALLEL -o $@
