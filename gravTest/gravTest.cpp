#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include "Vecmat.h"

/*	The following functions are heavily based on the works of Samual Pines:
 
	Pines, S., “Uniform Representation of the Gravitational Potential and its Derivatives,” AIAA Journal,
	Vol. 11, No. 11, Nov. 1973, pp. 1508–1511.
 
	And the examples found in:
	
	Normalization and Implementation of Three Gravitational Acceleration Models
	Randy A. Eckman
	Aaron J. Brown
	Daniel R. Adamo


*/


class TessGravProp
{
public:
	//TessGravProp();
	~TessGravProp();
	bool readGravModel(char* filename, int cutoff);
	Vector GetTessGrav(const Vector rposmax, const int maxDegree, const int maxOrder);
	void GenerateAssocLegendreMatrix(int maxDegree);

	static inline unsigned int NM(unsigned int n, unsigned int m) {return (n * n + n) / 2 + m;}
	static inline double KroneckerDelta(int m) { return (double)((m == 0) ? 1 : 0); };
	static inline double normLegendre(unsigned int n, unsigned int m) {return sqrt((2 - KroneckerDelta(m)) * (2 * (double)n + 1) * (std::tgamma(n - m + 1) / std::tgamma(n + m + 1)));}

	double refRad;
	double GM;
	unsigned int degree;
	unsigned int order;
	unsigned int normalized;
	double referenceLat;
	double referenceLon;
	double* C;
	double* S;
	double* R;
	double* I;
	unsigned long int numCoeff;

	double r, s, t, u;
	double rho;
	double* A;
	double a1, a2, a3, a4;
};


int main()
{
	TessGravProp gravityProperties;

	//char gravModelName[256] = "jgl075d1.sha";
	char gravModelName[256] = "jggrx_1500e_sha";


	unsigned int maxDegree = 0;

	std::cin >> maxDegree;
	

	Vector R = Vector(1738.0, 30.0, 0.0);

	bool isload = gravityProperties.readGravModel(gravModelName, maxDegree);
	

	auto start = std::chrono::high_resolution_clock::now();

	Vector A(0.0, 0.0, 0.0);

	if (isload){
		std::cout << "Loaded...Starting..." << std::endl;
		Vector A = gravityProperties.GetTessGrav(R, maxDegree, maxDegree);
	}

	
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
	std::cout << A << std::endl;
	std::cout << "Time Per Step:\t\t" << (double)duration.count() << std::endl;

	system("pause");
}

bool TessGravProp::readGravModel(char* filename, int cutoff)
{
	FILE* gravModelFile;
	char gravFileLine[512];
	bool isEOF = false;
	unsigned int linecount = 0;
	unsigned int maxLines = (cutoff * cutoff + cutoff) / 2 + cutoff;
	C = new double[1];
	S = new double[1];
	numCoeff = 0;

	C[0] = 0;
	S[0] = 0;

	std::cout.precision(16);

	fopen_s(&gravModelFile, filename, "rt");

	if (gravModelFile) {
		while (fgets(gravFileLine, 511, gravModelFile))
		{
			if (feof(gravModelFile))
			{
				break;
			}

			if (linecount == 0) {
				if (!sscanf(gravFileLine, " %lf , %lf , %*lf , %d , %d , %d , %lf , %lf \n",
					&refRad,
					&GM,
					&order,
					&degree,
					&normalized,
					&referenceLat,
					&referenceLon)) {
					return false;
				}
				numCoeff = linecount++;
			}
			else if (linecount <= maxLines) {
				double* tempC = new double[linecount];
				double* tempS = new double[linecount];

				memcpy(tempC, C, (linecount) * sizeof(double));
				memcpy(tempS, S, (linecount) * sizeof(double));

				delete[] C;
				delete[] S;

				C = tempC;
				S = tempS;


				if (!sscanf(gravFileLine, " %*d , %*d , %lf , %lf , %*lf , %*lf \n",
					&C[linecount - 1],
					&S[linecount - 1])) {
					return false;
				}
				numCoeff = linecount++;
			}
		}
		return true;
	}
	else {
		//some error goes into the log here
		return false;
	}
}

Vector TessGravProp::GetTessGrav(const Vector rpos, const int maxDegree, const int maxOrder)
{
	if ((maxDegree < 2) || (maxOrder < 2)) {
		Vector temp;
		temp.x = 0;
		temp.y = 0;
		temp.z = 0;
		return temp;
	}	

	r = rpos.length();
	s = rpos.x / r;
	t = rpos.y / r;
	u = rpos.z / r;

	rho = GM / (r*refRad);
	double rhop = refRad / r;

	R = new double[maxOrder+1];
	I = new double[maxOrder+1];

	R[0] = 0.0;
	I[0] = 0.0;
	R[1] = 1.0;
	I[1] = 0.0;

	for (int m = 2; m <= maxOrder; m++) {
		R[m] = s * R[m - 1] - t * I[m - 1];
		I[m] = s * I[m - 1] + t * R[m - 1];
	}

	double g1temp = 0.0;
	double g2temp = 0.0;
	double g3temp = 0.0;
	double g4temp = 0.0;

	double g1 = 0.0;
	double g2 = 0.0;
	double g3 = 0.0;
	double g4 = 0.0;

	int nmodel = 0;

	GenerateAssocLegendreMatrix(maxDegree);


	for (int n = 0; n <= maxDegree; n++) {
		

		double SM = 0.5;

		if (n > maxOrder)
			nmodel = maxOrder;
		else
			nmodel = n;



		for (int m = 0; m <= nmodel; m++) {
			double D = C[NM(n, m)] * R[m+1] + S[NM(n, m)] * I[m+1];
			double E = C[NM(n, m)] * R[m] + S[NM(n, m)] * I[m];
			double F = S[NM(n, m)] * R[m] - C[NM(n, m)] * I[m];

			double ALPHA = sqrt(SM*(n-m)*(n+m+1));

			g1temp = g1temp + A[NM(n, m)] * m * E;
			g2temp = g2temp + A[NM(n, m)] * m * F;
			g3temp = g3temp + ALPHA * A[NM(n, m)] * D;
			g4temp = g4temp + ((n + m + 1) * A[NM(n, m)] + ALPHA * u * A[NM(n, m+1)] * D);

		}
	}

	rho = rhop * rho;

	g1 = g1 + rho * g1temp;
	g2 = g2 + rho * g2temp;
	g3 = g3 + rho * g3temp;
	g4 = g4 + rho * g4temp;


	Vector gperturbed;

	gperturbed.x = (g1 - g4 * s);
	gperturbed.y = (g2 - g4 * t);
	gperturbed.z = (g3 - g4 * u);

	return gperturbed;
}

void TessGravProp::GenerateAssocLegendreMatrix(int maxDegree)
{
	A = new double[NM(maxDegree + 3, maxDegree + 3)]; //FIXME move to read function
	A[0] = sqrt(2.0);

	for (int m = 0; m <= (maxDegree+2); m++) {

		if (m != 0) {
			A[NM(m, m)] = sqrt(1. + (1. / (2. * (double)m))) * A[NM(m-1, m-1)]; // diagonal terms
		}

		if (m != (maxDegree + 1)) {
			A[NM(m + 1, m)] = sqrt(2. * (double)m + 3.) * u * A[NM(m, m)]; // off-diagonal terms
		}

		if (m < maxDegree) {
			for (int n = m + 2; n <= (maxDegree+2); n++) {
				double ALPHA_NUM = (2. * (double)n + 1.) * (2. * (double)n - 1.);
				double ALPHA_DEN = ((double)n - (double)m) * ((double)n + (double)m);
				double ALPHA = sqrt(ALPHA_NUM / ALPHA_DEN);
				double BETA_NUM = (2. * (double)n + 1.) * ((double)n - (double)m - 1.) * ((double)n + (double)m - 1.);
				double BETA_DEN = (2. * (double)n - 3.) * ((double)n + (double)m) * ((double)n - (double)m);
				double BETA = sqrt(BETA_NUM / BETA_DEN);


				A[NM(n, m)] = ALPHA * u * A[NM(n-1, m)] - BETA * A[NM(n-2, m)]; // remaining terms in the column
			}
		}
	}

	for (int n = 0; n <= (maxDegree+2); n++) {
		A[NM(n, 0)] = A[NM(n, 0)] * sqrt(0.5);
	}

	//for (int n = 0; n <= maxDegree+2; n++) {
	//	for (int m = 0; m <= n; m++) {
	//		std::cout << n << "\t" << m << "\t" << A[NM(n, m)] << std::endl;
	//	}
	//}
}

TessGravProp::~TessGravProp()
{
	delete[] C;
	delete[] S;
	delete[] A;
}