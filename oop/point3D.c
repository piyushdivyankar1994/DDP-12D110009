point3D* point3D_newPoint(float x, float y, float z){
	point3D* new = (point3D*)malloc(sizeof(point3D));
	new->x = x;
	new->y = y;
	new->z = z;
	return new;
}

point3D* point3D_origin(){
	point3D* new = (point3D*)malloc(sizeof(point3D));
	new->x = 0;
	new->y = 0;
	new->z = 0;
	return new;
}

point3D* point3D_addVectors(point3D* v1, point3D* v2){
	point3D* result = (point3D*)malloc(sizeof(point3D));
	result->x = v1->x + v2->x;
	result->y = v1->y + v2->y; 
	result->z = v1->z + v2->z;
	return result;
}

int point3D_isEqual(point3D* v1, point3D* v2) {
	if((v1->x - v2->x) < 1e-6 && (v1->y - v2->y) < 1e-6 && (v1->z - v2->z) < 1e-6) {
		return 1;
	}
	else return -1;
}

point3D* point3D_subtractVectors(point3D* v1, point3D* v2){
	point3D* result = (point3D*)malloc(sizeof(point3D));
	result->x = v1->x - v2->x;
	result->y = v1->y - v2->y; 
	result->z = v1->z - v2->z;
	return result;
}

void point3D_dispPoint(point3D* a){
	if(a == NULL) return;
	printf("(%.2f\t%.2f\t%.2f)\t\n", a->x, a->y, a->z);
	return;
}

point3D* point3D_indexToPoint3D_fcc(int index, parameter* p){
	point3D* returnPoint = point3D_origin();
	int motiff_position = index % p->atoms_per_site;          //type of atom
	
	returnPoint->x = (float)((index / p->atoms_per_site) % p->Nx);
	returnPoint->y = (float)((index / (p->atoms_per_site * p->Nx)) % p->Ny);
	returnPoint->z = (float)((index / (p->atoms_per_site * p->Nx * p->Ny)) % p->Nz);

	if(motiff_position == 0)                                        //if atom at first position(0,0,0) indicies are the global frame coordinates
	{
		return returnPoint;
	}
	else if (motiff_position == 1) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0, 0.5, 0.5));
	}
	
	else if (motiff_position == 2) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0.5, 0, 0.5));
	}
	
	else if (motiff_position == 3) {
		return point3D_addVectors(returnPoint, point3D_newPoint(0.5, 0.5, 0));
	}
	
	return returnPoint;
}

int point3D_point3DtoIndex(point3D* a, parameter* p){
	float fx = a->x - floor(a->x);
	float fy = a->y - floor(a->y);
	float fz = a->z - floor(a->z);
	
	int pos = ((int)floor(a->x) + p->Nx * (int)floor(a->y) + p->Nx * p->Ny * floor(a->z)) * 4;
	
	if(fx == 0 && fy == 0 && fz == 0){
		return pos;
	}
	else if(fx == 0)
		return pos + 1;
	else if(fy == 0)
		return pos + 2;
	else if(fz == 0)
		return pos + 3;
	return pos;
}

void point3D_periodicBoundaryTransform(point3D* k, parameter* p) {
	if(k->x < 0) {
		k->x = k->x + p->Nx;
	}
	else if(k->x > ((float)p->Nx - 0.5)){
		k->x = k->x - p->Nx;
	}
	if(k->y < 0) {
		k->y = k->y + p->Ny;
	}
	else if(k->y > ((float)p->Ny - 0.5)){
		k->y = k->y - p->Ny;
	}
	if(k->z < 0) {
		k->z = k->z + p->Nz;
	}
	else if(k->z > ((float)p->Nz - 0.5)){
		k->z = k->z - p->Nz;
	}
	//return k;
}

float point3D_magnitude(point3D* a){
	return  sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
}

float point3D_distAtoB(point3D* A, point3D* B){
	return point3D_magnitude(point3D_subtractVectors(A, B));
}
