double alnorm ( double x, bool upper );
double tfn ( double x, double fx );

extern "C" {  
	void psnc(double *out, double *x, double xi, double omega, double alpha, int n);
	void diffc(double *out, double *x, int n);
}

