/*
datos del CAMB
H=70
T=2.725
Ob=0.002
Om=0.275551
Ol=0.0.722449
HeliumFraction=0.24
MasslessNeutrinos=0
MassiveNeutrino=0
ComovingSoundSpeed=1
*/


// usando 1 nucleo
// 100^3, 4 escalas, 3 muestras = 34 minutos ~ 3 minutos por escala por muestra
// 200^3, 4 escalas, 1 muestras = 96:30 minutos ~ 23 minutos por escala por muestra
// 500^5, cada escala deberia demorar 6 horas (~3:40 horas usando 4 nucleos)

//	variables para el fftw


	fftw_complex *delta;						//	potencial graitacional
	fftw_complex *sx, *sy, *sz;			//	desplazamientos de Zeldovich, cond. iniciales
	fftw_plan dft,dfti;							//	plan principal
	fftw_plan dftx,dfty,dftz;				//	plan para las condiciones iniciales


