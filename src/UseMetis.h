#ifndef _MPI_
#pragma once
#include<mpi.h>
#include <cstddef> /* NULL */
#include <metis.h>
#include <iostream>
#include<fstream>
#include <vector>

class UseMetis
{
public:
	UseMetis();
	~UseMetis();
	void geteptreind();
	void partiGraph(int*,int*,int*,int,int);

	std::vector<idx_t> eptr;
	std::vector<idx_t> eind;
	idx_t nElements;// = eptr.size()-1;

	idx_t nn;

	idx_t objval;//存分区总通信量


};
#endif
