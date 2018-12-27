#include <iostream>
#include <fstream>
#include <cstdlib>
#include <Eigen/Dense>
#include "FMM2DTree.hpp"
#include <ctime>
#include <iomanip> 
#include <stdlib.h> 

int main(int argc, char* argv[]) {
	//int nLevels		=	atoi(argv[1]);
	int nChebNodes	=	atoi(argv[1]);
	double L			=	std::strtod (argv[2],NULL);//atoi(argv[2]);
	
	srand(time(NULL));
	double start, end;

	start	=	omp_get_wtime();

	pts2D dummy;
	dummy.x = 0.0;
	dummy.y = 0.0;
	

	user_Greens_function* mykernel		=	new user_Greens_function;
	mykernel->initialise(dummy,dummy,0.0,0.0,0.0,0.0,0.0);
	FMM2DTree<user_Greens_function>* A	=	new FMM2DTree<user_Greens_function>(mykernel, nChebNodes, L);

	A->set_Standard_Cheb_Nodes();
	A->get_interpolant_matrix();

	A->createTree();
	A->assign_list2_Neighbor_Interactions0();
	A->modifyTree();

	A->assign_list2_Neighbor_Interactions();

	A->assign_list1_list3_Interactions();
	

	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);


	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;


	start	=	omp_get_wtime();

	A->assemble_Operators_FMM();
	A->get_Transfer_Matrix();
		
	end		=	omp_get_wtime();
	double timeAssemble		=	(end-start);
	std::cout << std::endl << "Time taken to assemble the operators is: " << timeAssemble << std::endl;

	
	start	=	omp_get_wtime();

	A->evaluate_multipoles();

	end		=	omp_get_wtime();
	double timeAssignCharges=	(end-start);
	std::cout << std::endl << "Time taken to assemble the charges is: " << timeAssignCharges << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_All_M2M();
	end		=	omp_get_wtime();
	double timeM2M			=	(end-start);
	std::cout << std::endl << "Time taken for multipole to multipole is: " << timeM2M << std::endl;

	start	=	omp_get_wtime();

	A->evaluate_list2();
	A->evaluate_list3();
	A->evaluate_list4();

	end		=	omp_get_wtime();
	double timeM2L			=	(end-start);
	std::cout << std::endl << "Time taken for multipole to local is: " << timeM2L << std::endl;

	start	=	omp_get_wtime();
	A->evaluate_All_L2L();
	end		=	omp_get_wtime();
	double timeL2L			=	(end-start);
	std::cout << std::endl << "Time taken for local to local is: " << timeL2L << std::endl;

	start	=	omp_get_wtime();
	
	A->evaluate_list1();
	
	end		=	omp_get_wtime();
	double timeLeaf			=	(end-start);
	std::cout << std::endl << "Time taken for self and neighbors at leaf is: " << timeLeaf << std::endl;

	double totalTime	=	timeCreateTree+timeAssemble+timeAssignCharges+timeM2M+timeM2L+timeL2L+timeLeaf;
	
	double applyTime	=	timeM2M+timeM2L+timeL2L+timeLeaf;

	std::cout << std::endl << "Total time taken is: " << totalTime << std::endl;

	std::cout << std::endl << "Apply time taken is: " << applyTime << std::endl;

	std::cout << std::endl << "Total Speed in particles per second is: " << A->N/totalTime << std::endl;

	std::cout << std::endl << "Apply Speed in particles per second is: " << A->N/applyTime << std::endl;

	
	std::cout << std::endl << "Performing Error check..." << std::endl;

	/*srand(time(NULL));
	int nBox	=	rand()%A->nBoxesPerLevel[nLevels];
	std::cout << std::endl << "Box number is: j: 2; k: 1" << std::endl;
	int j=2, k=1;
	std::cout << std::endl << "Box center is: (" << A->tree[j][k].center.x << ", " << A->tree[j][k].center.y << ");" << std::endl;*/
	//std::cout << std::endl << "Error is: " << A->perform_Error_Check() << std::endl;

	start	=	omp_get_wtime();

	A->perform_Error_Check();
	
	end		=	omp_get_wtime();
	double errorTime	=	(end-start);

	std::cout << std::endl << "Time taken to compute error is: " << errorTime << std::endl;
	
	std::cout << std::endl;

}
