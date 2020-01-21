 

void crece(double &D_i,double &Dp_i,double &D_antes,double &Dp_antes){
cout<<"computing growth function"<<endl;
double y1,y2,a;//las variables del Runge Kutta
//FILE * C; 	C=fopen("crecimiento.dat","w+");
y1=0.00001;	y2=0.;	a=0.00001;

cout<<"calculando la funcion de crecimiento"<<endl;
FILE * FI;
FI=fopen("d.dat","w+");
int i=1;//contador para imprimir una muestra de la funcion de crecimiento
  while(a<=ai-0.5*aa){//ciclo en el factor de escala hasta llegar a (ai+da/2)
    RungeKutta(a,y1,y2,0.0000001); 
    a+=0.0000001;
fprintf(FI,"%lf %lf %lf \n",a,y1,y2*a*H(a));
//     if (fmod(i*1.,1000.)==0.){fprintf(C,"%le %le %le\n",a,y1,y2*a*H(a));}i++;
}
  Dp_antes=y2*a*H(a);
  D_antes=y1;

  while(a<=ai){//ciclo en el factor de escala hasta llegar a (ai)
    RungeKutta(a,y1,y2,0.0000005);
    a+=0.0000005;
//     if (fmod(i*1.,1000.)==0.){fprintf(C,"%le %le %le\n",a,y1,y2*a*H(a));}i++;
fprintf(FI,"%lf %lf %lf \n",a,y1,y2*a*H(a));
}
  Dp_i=y2*a*H(a);
  D_i=y1;

  while(a<=1.){//ciclo en el factor de escala hasta llegar a (1) para obtener la constante de normalizacion
    RungeKutta(a,y1,y2,0.00001);
    a+=0.00001;
fprintf(FI,"%lf %lf %lf \n",a,y1,y2*a*H(a));
//     if (fmod(i*1.,1000.)==0.){fprintf(C,"%le %le %le\n",a,y1,y2*a*H(a));}i++;
}
printf("\n\nnorma\n\n %lf\n",y1);
fclose(FI);

  Dp_i/=y1;				D_i/=y1;	//normalizando la funcion de crecimiento y su derivada
  Dp_antes/=y1;		D_antes/=y1;	//normalizando la funcion de crecimiento y su derivada

  cout<<"ai="<<a<<"  Di="<<D_i<<"   Dpi="<<Dp_i<<""<<endl;//crecimiento y derivada para las condiciones iniciales
  cout<<"ai-="<<a<<"  D(i-d)="<<D_antes<<"   Dp(i-d)="<<Dp_antes<<""<<endl<<endl;
  cout<<"delta a="<<aa<<endl;
  cout<<"ai="<<ai<<endl;
//fclose(C);


}
