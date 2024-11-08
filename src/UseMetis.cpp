#ifndef _MPI_

#include "UseMetis.h"

void UseMetis::geteptreind()
{
	//读入分块文件
	std::ifstream data_eptr("eptr.txt");
	if (!data_eptr.is_open())
	{
		std::cout << "Error File Path eptr!!!" << std::endl;
		// system("pause");
	}
	double point_eptr;
	data_eptr >> nn;//总节点数
	while (data_eptr >> point_eptr)
		eptr.push_back(point_eptr);

	data_eptr.close();


	std::ifstream data_eind("eind.txt");
	if (!data_eind.is_open())
	{
		std::cout << "Error File Path eind!!!" << std::endl;
		// system("pause");
	}
	double point_eind;
	while (data_eind >> point_eind)
		eind.push_back(point_eind);

	data_eind.close();

	nElements = eptr.size() - 1;

}
void UseMetis::partiGraph(int *pE2d,int *pE3d,int*pV,int nL,int nParts)
{
	//partElements与partNodes中的nParts从0开始
	int ret = METIS_PartMeshNodal(&nElements, &nn, eptr.data(), eind.data(), NULL, NULL, &nParts, NULL, NULL, \
		&objval, pE2d, pV);

	//std::cout << "Metis is OK! " << std::endl;
	//std::cout << "objval: " << objval << std::endl;

	//DG只需要存每个单元在哪个分区即可，对应K2d、K3d的一个全局变量
	for (int part_i = 0; part_i < nElements; part_i++) {
		//std::cout << "The element " << part_i + 1 << " is partied in nPart " << pE2d[part_i] << std::endl;
		for (int part_j = 0; part_j < nL; part_j++) {
			pE3d[part_i * nL + part_j] = pE2d[part_i];
		}

	}

}

//void UseMetis::sendrecmpi(double *fphys, int *Np, int*K, int Nfield)
//{
//	MPI_Status status;
//	int count = (*Np)*(*K)*Nfield;
//	int size;
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	double buff[(*Np)*(*K)*(Nfield)];
//	if (id != 0)
//	{
//		int send = MPI_Send(fphys, count, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
//	}
//	//std::cout <<"here"<<id<< std::endl;
//	if (id == 0)
//	{
//		for (int i = 1; i < (size); i++)
//		{
//			int rec = MPI_Recv(buff, count, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
//			for (int k = 0; k < *K; k++)
//			{
//				if (part_r[k] != i) continue;
//				for (int n = 0; n < *Np; n++)
//				{
//					int sk = k * (*Np) + n;
//					fphys[sk] = buff[sk];
//				}
//			}
//		}
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	if (id == 0)
//	{
//		for (int i = 1; i < (size); i++)
//		{
//			int rec = MPI_Send(fphys, count, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
//
//		}
//	}
//	if (id != 0)
//	{
//		int rec = MPI_Recv(fphys, count, MPI_DOUBLE, 0, id, MPI_COMM_WORLD, &status);
//
//	}
//
//}

UseMetis::UseMetis()
{
	geteptreind();

}


UseMetis::~UseMetis()
{

}

#endif