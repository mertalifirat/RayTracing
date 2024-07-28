all:
	g++ *.cpp -O3 -pthread -o raytracer -std=c++11
debug:
	g++ *.cpp -g -O3 -pthread -o raytracer -std=c++11
