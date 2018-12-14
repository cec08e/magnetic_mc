default: heisenberg3d_nlayer.c plot_h3d_nlayer.py
	gcc -fPIC -shared -o heisenberg3d_nlayer.so -g -I /usr/include -L /usr/include heisenberg3d_nlayer.c -lm -lgsl -lgslcblas 
	python3 plot_h3d_nlayer.py
