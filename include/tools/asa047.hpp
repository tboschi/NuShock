typedef double (*function_t)( double x[] , void *userdata, int sign);

void nelmin ( function_t fn, int n, void *userdata, int sign, double start[], double xmin[], 
	      double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
	      int *icount, int *numres, int *ifault );

void timestamp ( );
