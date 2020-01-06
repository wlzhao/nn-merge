# On the Merge of k-NN Graph

## About the codes

This is the implementation of three k-NN graph merge algorithms, S-Merge,J-Merge and H-Merge. The description about the algorithms could be found in "**On the Merge of k-NN Graph**".

* smerge.h smerge.cpp		  Symmetric Merge  
* jmerge.h jmerge.cpp		  Joint Merge
* lynndescent.h lynndescent.cpp	  Hierarchical k -NN graph Construction via J-Merge
* lysearch.h lysearch.cpp         NN search task based on H-Merge

## Compilation

* 1. The codes should be compiled with g++ 5.4.0 or later. 
* 2. Add following requirement in your linker options:
	* -lboost_timer
	* -lboost_system
* 3. In addition, boost library 1.58 is required. For higher version of boost may not be compatible. 


To compile it, user only needs to run "make Makefile" under the source code directory

## Run
To run the above mentioned methods, user is recommended to call "test()" function in "main.cpp". In the "test()" function, four interfaces, namely S-Merge, J-Merge, H-Merge and NN search algorithm, are provided. Furthermore, Dataset SIFT100K is attached with this package


Author: Peng-Cheng Lin, Wan-Lei Zhao
Date: Jan.-6-2020
