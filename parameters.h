const int nt=40;			//	number of threads for paralelization
// using 4 threads
// 100^3 particles and cells, 4 boxes, 3 statistical samples = 34 minutes ~ 3 minutes per escale per samples
// 200^3 particles and cells, 4 boxes, 1 statistical samples = 96:30 minutes ~ 23 minutes por escale por sample
// 128^3 particles and 128^3 cells ~8 minutes (just 1 box and 1 sample)
// 128^3 particles and 256^3 cells ~16 minutes (just 1 box and 1 sample)

///////////////////////////			Z inicial, final, paso...
const double ai=0.02;					//	initial scale factor
const double af=1.0;					//	final scale factor
const double aa=0.01;					//	steep size in a


///////////////////////////			parámetros cosmológicos
const double Ol=0.68887;				//	fracción de energia oscura	//0.6825;
const double Om=0.311051;				//	fracción de materia		//0.3175;
const double Ok=0.;				//	fracción de curvatura hoy
const double Og=0.;				//	fracción de radiacion				//0.0000475


char file_CAMB[]="planck_2018_matterpower_49.dat";
char file_IC[]="initial_conditions_from_2LPT.dat";

////////////////////////////		numero de particulas, celdas, y derivados
int Lbox=512;				//	box side in Mpc/h
int n=256;				//	number of cells
int np=256;				//	numer of particles to the power -3

int GADGET_IC = 0;			//	Initial conditions used by GADGET
