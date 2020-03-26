
INC_DIR := -I/usr/local/opt/petsc/include
LIB_DIR := -L/usr/local/opt/petsc/lib
#LIB_DIR := -L/home/alain/Bureau/Projet_2660/lib_petsc/lib -Wl,-rpath=/home/alain/Bureau/Projet_2660/lib_petsc/lib
LIB := -lpetsc.3.12

CXX_FLAGS := -O0 -Ofast -lstdc++ #-std=c++14 #-g

#Compilation 
all : 
	gcc -o project -x c poisson.c -x c++ project.c airfoil.cpp -lm -lmpi $(CXX_FLAGS) $(LIB_DIR) $(LIB) $(INC_DIR)

#Delete of executable file
clean :
	rm projet

#Delete of results
clean_txt :
	rm -vf results/*.txt results/P-* results/U-* results/V-* results/Reh-* results/Vtx-* results/Rehw-* results/Div-*
