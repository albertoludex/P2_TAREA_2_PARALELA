#include <stdlib.h>
#include <stdio.h>
#include<string.h>
#include <math.h>
#include<time.h>
#include <omp.h>
#define TAMANO 15000

double** matrix_create(int ancho, int largo){
    int i;
    double** matriz = (double**) malloc (largo * sizeof(double*)); //Se declara un array de punteros del tamaño de largo
    for(i=0;i<largo;i++){
        matriz[i]= (double*) malloc (ancho * sizeof(double)); //Por cada puntero del array anterior se reserva memoria para una array de unsigned chars de tamaño ancho
    }
    
    return matriz; 
    
}
void matrixxvector(double** matriz, double* Vorg, double* Vdes, int N,int largo, int div, int iam){  

    int i,j;
    double var=0;
    for(i=iam*div;i<largo;i++){ //Fila
        var=0;
        for(j=0;j<N;j++){//columna 
        var+=Vorg[j]*matriz[i][j];
        }
        Vdes[i]=var; //utilizamos el mismo iterador que utilizamos para recorrer la matriz fila a fila para recorrer el vector destino elemento a elemeto
    }
}
//Funcion que crea vectores de tamaño N y los inicializa a un valor ini que le pasaremos por valor
double* vector_create(int N, double ini){
    double* vector;
    int i;
    vector= (double*) malloc (N * sizeof(double)); 
    for(i=0;i<N;i++)
    vector[i]=ini; 

    return vector;
}

void matrix_rellena(double** matriz, int N, char* fichero){
    int i,j;
    time_t t;
    FILE *file;
    if (file = fopen (fichero, "rb")){
        for(i=0;i<N;i++){        
            fread(matriz[i],sizeof(double), N, file);           
        }
        fclose(file);
        
    }
    else{
        //CODIGO DE COMPROBACION
        
        for (i=0;i<N;i++){

            for (j=0;j<N;j++){
                if (i == j){
                    matriz[i][j] = 1.0;
                }
                else if(i > j){
                    matriz[i][j]=(double)-50*(i+1)*(j+1)/(N+N); //Triangular superior
                }
                else{
                    matriz[i][j]=(double)50*(i+1)*(j+1)/(N+N); //Triangular inferior
                }
                

            }
        
        
        }
        
        /*
    
        for(i=0;i<N;i++){ //Fila matriz
            for(j=0;j<N;j++){ //recorre filas
                if(i==j){ //diagonal
                    matriz[i][j]=1;
                }
                else if(i>j){//Diagonal inferior
                    matriz[i][j]= rand() / ((double) RAND_MAX)*-50;
                }
                else{//Diagonal superior
                    matriz[i][j]= rand() / ((double) RAND_MAX)*50; 
                }
            }
        }*/
        
        file = fopen (fichero, "wb"); 
        for(i=0;i<N;i++){        
            fwrite(matriz[i],sizeof(double), N, file);
        }
        
    }
}

void dividiendo(double* vector,int N,double divisor,int iam,int div){
    int i;
    for(i=iam*div;i<N;i++){
        vector[i]=vector[i]/divisor;
    }
}

double calcular_abs(double* vect, int N){
int i;
    double temp=0;
    double valor_literal=25.0;
    

    //Recorre el vector y encuentra el máximo valor absoluto
    for(i=0;i<N;i++){
        if(sqrt(pow(vect[i],2))>temp){
            temp=sqrt(pow(vect[i],2));
        }
    }
    //divide el vector por el maximo valor absoluto
    for(i=0;i<N;i++){

        //lo comparamos con el valor literal 
        if (vect[i]>valor_literal)
        {
        vect[i]/temp;
        }
        
    }
    return temp;

	

}


int comp_fichero (char *fichero){
    FILE *f= fopen(fichero,"r");
    if (f==NULL){
        return 0;
        //Si no es un fichero devuelve 0 
    }
    fclose(f);
    return 1;
        //Si es un fichero devuelve 1
}

void ImpResultados(int itera,int nproces,double* resultados, double tt,char* fichero){
    char nombre[50];
    int contador=0;
    FILE *f;
    int i;
    
    sprintf(nombre,"Resultados_itera_%d.txt",itera);
     while(comp_fichero(nombre)==1){
        contador++;
        sprintf(nombre,"Resultados_itera_%d(%d).txt",itera,contador);
    }
   

    f=fopen(nombre,"w");
    fprintf(f,"Numero de iteraciones: %d\n\n",itera);
    fprintf(f,"Numero de hilos: %d\n\n",nproces);
    fprintf(f,"La primera iteracion no genera valor absoluto\n");
    for(i=0;i<itera-1;i++){
        fprintf(f,"Mayor absoluto iteración %d: ",i+2);
        fprintf(f,"%f  ",resultados[i]);
        fprintf(f,"\n");
    }

    fprintf(f,"\n\nTiempo paralelo total es %f\n",tt);
    fprintf(f,"El fichero de entrada y de salida es %s\n",fichero);
    

    fclose(f);
}  

int main(int argc, char* argv[]){
    srand(time(NULL));
    int i, j ,k,iam;
    int hilos = atoi(argv[3]);
    double absolut = 0;
    double var;
    int itera=atoi(argv[1]);
    double* absolutos= vector_create(hilos,0);
    double* resultados= vector_create(itera-1,0);
    double** M;
    double* V0 = vector_create(TAMANO,1);
    double* V1 = vector_create(TAMANO,1);
    char* fichero=argv[2];
    double t0, t1, t2, t3, t4;
    double tp=0;

    M= matrix_create(TAMANO,TAMANO);
    matrix_rellena(M,TAMANO,fichero);
    t0=omp_get_wtime();
    #pragma omp parallel num_threads(hilos) shared(M,hilos,absolutos,V1,V0) private(iam,var,j,i)
    {

        #pragma omp for schedule(static,(TAMANO/hilos))
        #pragma omp fo schedule (static, 2) reduction(+:var)

        for(i=0;i<TAMANO;i++){ //Fila
            var=0;
            for(j=0;j<TAMANO;j++){//columna 
                var+=V0[j]*M[i][j];
            }
            V1[i]=var;//utilizamos el mismo iterador que utilizamos para recorrer la matriz fila a fila para recorrer el vector destino elemento a elemeto
            }
    
    }
    
    t1=omp_get_wtime();
    for(k=1;k<itera;k++){
        if((k%2)==1){
            t3=omp_get_wtime();
            #pragma omp parallel num_threads(hilos) shared(M,hilos,absolut,absolutos,V0,V1,resultados,k) private(iam,var,i,j)
            {
        
                iam=omp_get_thread_num();
                #pragma omp for schedule(static,(TAMANO/hilos))
                        #pragma omp fo schedule (static, 2) reduction(+:var)

                    for(i=0;i<TAMANO;i++){ //Fila
                        var=0;
                        for(j=0;j<TAMANO;j++){//columna 
                            var+=V1[j]*M[i][j];
                        }
                    V0[i]=var; //utilizamos el mismo iterador que utilizamos para recorrer la matriz fila a fila para recorrer el vector destino elemento a elemeto
                    }
                #pragma omp for schedule(static,(TAMANO/hilos))
                        #pragma omp fo schedule (static, 2) reduction(+:absolutos)

                    for(i=0;i<TAMANO;i++){
                        if(sqrt(pow(V0[i],2))>absolutos[iam]){
                        absolutos[iam]=sqrt(pow(V0[i],2));
            
                        }
        
                    }  
                   
                
                #pragma omp barrier
                #pragma omp master
                {
                    
                    
                    for(i=0;i<hilos;i++){
                        if(sqrt(pow(absolutos[i],2))>absolut){
                        absolut=sqrt(pow(absolutos[i],2));
                        }
                    }
                    resultados[k-1]=absolut; 
                    

                }
                 #pragma omp barrier 
                 
                #pragma omp for schedule(static,(TAMANO/hilos))
                for(i=0;i<TAMANO;i++){
                    V0[i]=V0[i]/absolut;
                }

            }

            t4=omp_get_wtime();
            tp+=t4-t3;
            for(i=0;i<hilos;i++){
                        absolutos[i]=0;
                    }
            absolut=0;
        }
        else{
            t3=omp_get_wtime();
            
            #pragma omp parallel num_threads(hilos) shared(M,hilos,absolut,absolutos,V0,V1,resultados,k) private(iam,var,i,j)
            {
        
                iam=omp_get_thread_num();
                #pragma omp for schedule(static,(TAMANO/hilos))
                        #pragma omp fo schedule (static, 2) reduction(+:var)

                    for(i=0;i<TAMANO;i++){ //Fila
                        var=0;
                        for(j=0;j<TAMANO;j++){//columna 
                            var+=V0[j]*M[i][j];
                        }
                    V1[i]=var; //utilizamos el mismo iterador que utilizamos para recorrer la matriz fila a fila para recorrer el vector destino elemento a elemeto
                    }
                
                #pragma omp for schedule(static,(TAMANO/hilos))
                    for(i=0;i<TAMANO;i++){
                        if(sqrt(pow(V1[i],2))>absolutos[iam]){
                        absolutos[iam]=sqrt(pow(V1[i],2));
            
                        }
        
                    }  
                 
                 
                #pragma omp barrier
                #pragma omp master
                {
                    
                    
                    
                    
                     for(i=0;i<hilos;i++){
                        if(sqrt(pow(absolutos[i],2))>absolut){
                        absolut=sqrt(pow(absolutos[i],2));
                        }
                    }
                    
                   resultados[k-1]=absolut; 
                    

                }
                 #pragma omp barrier
                 

                #pragma omp for schedule(static,(TAMANO/hilos))
                for(i=0;i<TAMANO;i++){
                    V1[i]=V1[i]/absolut;
                }
            }

            t4=omp_get_wtime(); 
            tp+=t4-t3;
            for(i=0;i<hilos;i++){
                        absolutos[i]=0;
                    }
            absolut=0;
        }   
        
    }

    ImpResultados(itera,hilos,resultados,(t1-t0)+tp,fichero);
    free(V0);
    free(M);
    free(V1);
    free(absolutos);
    free(resultados);
}