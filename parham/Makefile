default : libparham.so

lib%.so : %.cc
	g++ -fopenmp $< -shared -fPIC -o $@
