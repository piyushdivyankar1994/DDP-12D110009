#include "test.h"

void test_new_parameters(){
  parameter* new = NULL;
  new = new_parameters("parameter1.txt");
  print_parameters(new);
  free(new);
}

void test_dataRetrival() {
  double r = 4.123;
  binEAMpot* data = NULL;
  rdf*    nn = (rdf*)malloc(sizeof(rdf));
  eDen_df* a = (eDen_df*)malloc(sizeof(eDen_df));

  eam_data_read(&data, "file_list.txt", "Al", "Ni");
  nn = rdf_radius_retrive(data, r);
  printf("%d\n", nn->index);
  printf("%le\n", nn->p11);
  printf("%le\n", nn->p12);
  printf("%le\n", nn->p22);
  printf("%le\n", nn->eDen1);
  printf("%le\n", nn->eDen2);

  a = eDen_df_charge_density_retrive(data, 0.01);
  printf("\n");
  printf("%le\n", a->eDen);
  printf("%le\n", a->embed1);
  printf("%le\n", a->embed2);
}

void test_eam_data_read(){
  binEAMpot* data = NULL;
  eam_data_read(&data, "file_list.txt", "Al", "Ni");

  printf("%d\n", data->no_of_files);
  printf("%s\n", data->atom1);
  printf("%s\n", data->atom2);
}

void test_point3D_indexToPoint3D_fcc(){
	parameter* p = _defaultFCCparameter();
	printf("Generating 10 random numbers from 0 to 3999 and displaying corresponding points\n");
	for(int i = 0; i < 15; i++){
		int index = rand() % p->no_of_atoms;
		printf("%d\t=\t", index);
		point3D_dispPoint(point3D_indexToPoint3D_fcc(index, p));
	}
}

void test_point3D(){
	point3D* a1 = point3D_origin();
	point3D_dispPoint(a1);
	a1 = point3D_newPoint(1, 1.5, 4.3);
	point3D_dispPoint(a1);
	point3D* a2 = point3D_newPoint(1.1, 3.4, 5.3);
	point3D_dispPoint(a2);
	point3D* sum = point3D_addVectors(a1, a2);
	point3D_dispPoint(sum);
	point3D* diff = point3D_subtractVectors(a1, a2);
	point3D_dispPoint(diff);
	printf("Distance Test from point3D_origin\n");
	a1 = point3D_newPoint(1, 1, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(1, 1.5, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(1.4, 1, 1);
	printf("%f\n", point3D_magnitude(a1));
	a1 = point3D_newPoint(0, 0, 5);
	printf("%f\n", point3D_magnitude(a1));

	printf("Distance between two points\n");
	for(int i = 0; i < 3; i++){
		a1 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
		a2 = point3D_newPoint(rand() / (float) RAND_MAX, rand() / (float) RAND_MAX, rand() / (float) RAND_MAX);
		point3D_dispPoint(a1);
		point3D_dispPoint(a2);
		printf("Distance = %f\n\n", point3D_distAtoB(a1, a2));

	}

	return;
}

void test_defaultFCCparameter(){
	parameter* new = _defaultFCCparameter();
	print_parameters(new);
}


void test_Sn_fcc_readNeighbours_fromFile(){
	Sn_fcc* new = Sn_fcc_readNeighbours_fromFile("file_list_neighbours.txt");
	//print_Neighbours(new);
	printf("Everything is continuoulsy stored so if i goes from 0 \
	to sum(indicies) then s1n to s7n can be accessed as s1n[i] \
	this is demonstrated by printing 13 points\n ");
	for(int i = 0; i < 13; i++){
		point3D_dispPoint(&(new->s1n[i]));
	}
}


void test_AtomicMatrixRead(){
	parameter* p = _defaultFCCparameter();
	int *a = atomicMatrixRead("out.txt", p);
	print_AtomicMatrix(a, 20, 50);
}


void test_energyAtIndexFCC() {
	int index = 100;
	binEAMpot* data = NULL;
	eam_data_read(&data, "file_list.txt", "Al", "Ni");
	parameter* p = _defaultFCCparameter();
	int *a = atomicMatrixRead("out.txt", p);
	//test_AtomicMatrixRead();
	Sn_fcc* new = _defaultFCCNeighbours();
	double e = energyAtIndexFCC(index, a, data, p, new);
	printf("at index = %d energy = %f\n", index, e);

}

void test_point3D_point3DtoIndex(){
	parameter* p = _defaultFCCparameter();
	point3D* p1 = point3D_newPoint(2, 0, 0);
	point3D* p2 = point3D_newPoint(0.5, 0.5, 1);
	point3D* p3 = point3D_newPoint(2.0, 0.5, 0.5);
	point3D* p4 = point3D_newPoint(0.5, 1, 0.5);

	printf("%d\t", point3D_point3DtoIndex(p1, p));
	point3D_dispPoint(p1);
	printf("%d\t", point3D_point3DtoIndex(p2, p));
	point3D_dispPoint(p2);
	printf("%d\t", point3D_point3DtoIndex(p3, p));
	point3D_dispPoint(p3);
	printf("%d\t", point3D_point3DtoIndex(p4, p));
	point3D_dispPoint(p4);
}

void test_point3D_periodicBoundaryTransform(){
	point3D* new = point3D_newPoint(-10.5, 11.5, 16.5);
	parameter* p = _defaultFCCparameter();
	printf("test function for periodic boundary transform\n");
	point3D_dispPoint(new);
	//point3D* test =
	point3D_periodicBoundaryTransform(new, p);
	point3D_dispPoint(new);
}

void test_energyInMatrix() {
  binEAMpot* data = NULL;
  eam_data_read(&data, "file_list.txt", "Al", "Ni");
  parameter* p = _defaultFCCparameter();
  int *a = atomicMatrixRead("out.txt", p);
  //test_AtomicMatrixRead();
  Sn_fcc* new = _defaultFCCNeighbours();

  /// Example:
  double *energyMap = NULL;
  energyInMatrix(&energyMap, a, data, p, new);
  printEnergyMap(energyMap, 100, 200);
  free(energyMap);

}

void test_deltaEnergyMatrix() {
  binEAMpot* data = NULL;
  eam_data_read(&data, "file_list.txt", "Al", "Ni");
  parameter* p = _defaultFCCparameter();
  int *a = atomicMatrixRead("out.txt", p);
  //test_AtomicMatrixRead();
  Sn_fcc* new = _defaultFCCNeighbours();

  /// Example:
  double *energyMap = NULL;
  deltaEnergyMatrix(&energyMap, a, data, p, new);
  printEnergyMap(energyMap, 100, 200);
  free(energyMap);
}
