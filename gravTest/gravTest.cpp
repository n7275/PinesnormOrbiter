#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include "Vecmat.h"

class TessGravProp
{
public:
	//TessGravProp();
	~TessGravProp();
	bool readGravModel(char* filename, int cutoff);
	Vector GetTessGrav(const Vector rposmax, const double maxDegree);
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
	unsigned long int numCoeff;

	double r, s, t, u;
	double rho, R, I, A00;
	double* A;
	double a1, a2, a3, a4;
};


int main()
{
	TessGravProp gravityProperties;

	//char gravModelName[256] = "jgl075d1.sha";
	char gravModelName[256] = "jggrx_1500e_sha";

	//if(readGravModel(gravityProperties, gravModelName, 300))
	//{
	//	for (unsigned int i = 0; i < gravityProperties.numCoeff; i++)
	//	{
	//		std::cout << i << "   " << gravityProperties.C[i] << "   " << gravityProperties.S[i] << std::endl;
	//	}
	//}

	unsigned int maxDegree = 0;

	std::cin >> maxDegree;
	

	/*Vector R = Vector(1738.0, 0.0, 0.0);
	Vector A = gravityProperties.GetTessGrav(R, maxDegree);*/

	gravityProperties.u = 1;

	//bool isload = gravityProperties.readGravModel(gravModelName, maxDegree);
	gravityProperties.GenerateAssocLegendreMatrix(maxDegree);
	int n, m;
	std::cin >> n;
	std::cin >> m;;

	std::cout << gravityProperties.A[gravityProperties.NM(n, m)];

	//auto start = std::chrono::high_resolution_clock::now();

	
	//auto stop = std::chrono::high_resolution_clock::now();
	//auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	//std::cout << "Time Per Step:\t\t"(double)duration.count() << std::endl;

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
		GenerateAssocLegendreMatrix(cutoff);
		return true;
	}
	else {
		//some error goes into the log here
		return false;
	}
}

Vector TessGravProp::GetTessGrav(const Vector rpos, const double maxDegree)
{
	r = rpos.length();
	s = rpos.x / r;
	t = rpos.y / r;
	u = rpos.z / r;

	rho = GM / r;
	R = 1;
	I = 0;


	for (int n = 0; n <= maxDegree; n++) {
		for (int m = 0; m <= n; m++) {
			rho = pow(refRad / r, n + 1);
			
			a1 = 0;
			a2 = 0;

			R = (s * R) - (t * I);
			I = (s * I) + (t * R);

			a3 = 0;
			a4 = 0;
		}
	}

	return Vector();
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

	for (int n = 0; n <= maxDegree+2; n++) {
		for (int m = 0; m <= n; m++) {
			std::cout << n << "\t" << m << "\t" << A[NM(n, m)] << std::endl;
		}
	}
}

TessGravProp::~TessGravProp()
{
	delete[] C;
	delete[] S;
	delete[] A;
}