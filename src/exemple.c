#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include "math.h"

int Bas[3][3] = {
                  {1,1,1},
                  {1,1,1},
                  {1,1,1}
};

int xGradient[3][3] = {
                  {-1,0,1},
                  {-2,0,2},
                  {-1,0,1}
};

int yGradient[3][3] = {
                  {-1,-2,-1},
                  {0,0,0},
                  {1,2,1}
};

int Gaussian_filter[3][3] = {
                  {1,2,1},
                  {2,4,2},
                  {1,2,1}
};

int xCanny[3][3] = {
                  {-1,0,1},
                  {-1,0,1},
                  {-1,0,1}
};

int yCanny[3][3] = {
                  {-1,-1,-1},
                  {0,0,0},
                  {1,1,1}
};


int SEUIL = 30;
int NB_IMAGE = 967;

byte **dilatation(byte **m,long nrl,long nrh,long ncl,long nch);
byte **erosion(byte **m,long nrl,long nrh,long ncl,long nch);
byte **RGB_B(rgb8 **m,long nrl,long nrh,long ncl,long nch);
byte **ouverture(byte **m,long nrl,long nrh,long ncl,long nch,int iteration);
byte **binarisation(byte **m,long nrl,long nrh,long ncl,long nch,int seuil);
byte **fermeture(byte **m,long nrl,long nrh,long ncl,long nch,int iteration);

//Fonction initialisant dynamiquement une matrice
int **InitTab(int taillex,int tailley){

    int i;
    int **matrice = NULL;

    //allocation d'un tableau de tableaux d'entiers
    matrice = malloc ( taillex * sizeof(int *) );    

    for ( i = 0 ; i < taillex ; i ++ )
    {
        //allocation d'un tableau de tableau 
        matrice[i] = malloc ( tailley * sizeof(int *) ); 
    }
    return matrice;
}

/***  EXERCICE TP2 ***/

    




/***  EXERCICE TP1 ***/
//Fonction representant la methode de Harris
void Harris(char *path,char *imagename,char *pathfinal,char *pathFirstImage,int lambda){
    byte **matrice_temp;
    byte **matrice_difference;
    byte **matrice_gaussian;
    byte **matrice_harris;
    rgb8 **I;
    long nrh,nrl,nch,ncl;


    I=LoadPPM_rgb8matrix(pathFirstImage,&nrl,&nrh,&ncl,&nch);
    
    matrice_difference = bmatrix(nrl,nrh,ncl,nch);
    matrice_difference = RGB_B(I,nrl,nrh,ncl,nch);    

	int **tab = NULL;
	tab = InitTab(nrh*nch,256);

    rgb8 **R;
    int i,j,k;
    char str[50] = "\0";
    char num[4];
    char dest[50] = "\0";

	// On initialise le tableau à 0
	for(i = 0;i<nrh*nch;i++){
		for(j = 0; j < 256;j++){
			tab[i][j] = 0;
			//printf("i = %d && j = %d \n",i,j);
		}
	}
	printf("INITIALISATION DONE \n");

    for(i = 1;i < NB_IMAGE; i++){
        // On obtient le bon nom du fichier pour i :
        memset (str, 0, sizeof (str));
        strcat(str,path);
        strcat(str,imagename);
        if(i < 10){
            strcat(str,"00");
        }
        else if((i < 100)&&(i>=10)){
            strcat(str,"0");
        }

        snprintf(num, 4, "%d", i);

        strcat(str,num);
        strcat(str,".ppm");

        // On fait l'algorithme
        
        R=LoadPPM_rgb8matrix(str,&nrl,&nrh,&ncl,&nch);
        matrice_temp = bmatrix(nrl,nrh,ncl,nch);
        matrice_temp = RGB_B(R,nrl,nrh,ncl,nch);

        matrice_gaussian = bmatrix(nrl,nrh,ncl,nch);
        matrice_harris = bmatrix(nrl,nrh,ncl,nch);

        int GradientX = 0;
        int GradientY = 0;
        int Gaussian = 0;

        int compteur_c,compteur_l,ligne,colonne;
        //On fait notre filtre gaussien
        for(j = nrl;j < nrh;j++){
		    for(k = ncl;k < nch;k++){
				if((j == nrl)||(j == nrh - 1) || (k == ncl) || (k == nch - 1)){
                    matrice_gaussian[j][k] = 128;
                }
                else{
                    compteur_c = 0;
                    compteur_l = 0;

                    for(ligne = j-1;ligne<=j+1;ligne++){
                        for(colonne = k-1;colonne<=k+1;colonne++){
                            Gaussian = Gaussian + matrice_temp[ligne][colonne] * Gaussian_filter[compteur_l][compteur_c];
     
                            //On incrémente les compteurs
                            if(compteur_c ==2){
                                compteur_c = 0;
                                compteur_l++;
                            }  
                            else{
                                compteur_c++;
                            }
                        }
                    }
                    Gaussian = Gaussian/16;
                    matrice_gaussian[j][k] = Gaussian;
                    Gaussian = 0;
                }
            }
        }
        printf("Gaussian done\n");  
        
        int temp;
        //Calcul de gradient à partir de la matrice_gaussian
        for(j = nrl;j < nrh;j++){
            for(k = ncl;k < nch;k++){
                if((j == nrl)||(j == nrh - 1) || (k == ncl) || (k == nch - 1)){
                    matrice_harris[j][k] = 128;
                }
                else{
                    compteur_c = 0;
                    compteur_l = 0;

                    for(ligne = j-1;ligne<=j+1;ligne++){
                        for(colonne = k-1;colonne<=k+1;colonne++){
                            GradientX = GradientX + matrice_gaussian[ligne][colonne] * xCanny[compteur_l][compteur_c];
                            GradientY = GradientY + matrice_gaussian[ligne][colonne] * yCanny[compteur_l][compteur_c];
     
                            //On incrémente les compteurs
                            if(compteur_c ==2){
                                compteur_c = 0;
                                compteur_l++;
                            }  
                            else{
                                compteur_c++;
                            }
                        }
                    }

                    temp = (GradientX*GradientX)*(GradientY*GradientY) - GradientX*GradientY - lambda*(((GradientX*GradientX)+(GradientY*GradientY))*((GradientX*GradientX)+(GradientY*GradientY)));
                    if(temp >255){
                        matrice_harris[j][k] = 255;
                    }
                    else{
                        matrice_harris[j][k] = matrice_temp[j][k];
                    }
                    GradientX = 0;
                    GradientY = 0;
                }
            }
        }
        printf("Harris done\n");  
   

    //On envoie dans un dossier
    memset (dest, 0, sizeof (dest));
    strcat(dest,pathfinal);
    strcat(dest,"Harris");
    strcat(dest,num);
    strcat(dest,".pgm");    
    printf("dossier = %s\n",dest);
    SavePGM_bmatrix(matrice_harris,nrl,nrh,ncl,nch,dest);
    }
}

//Fonction permettant de faire l'étiquettage ////A faire\\\\



/**
* Fonction permettant de faire la différence par rapport à une image de fond
*/
void Exercice6(char *path,char *imagename,char *pathfinal,char *pathFirstImage){
    byte **matrice_temp;
    byte **matrice_fond;
    byte **matrice_difference;
    

    rgb8 **R;
    
    int i,j,k;
    long nrh,nrl,nch,ncl;


    char str[50] = "\0";
    char num[4];
    char dest[50] = "\0";
    for(i = 1;i < NB_IMAGE; i++){
        // On obtient le bon nom du fichier pour i :
        memset (str, 0, sizeof (str));
        strcat(str,path);
        strcat(str,imagename);
        if(i < 10){
            strcat(str,"00");
        }
        else if((i < 100)&&(i>=10)){
            strcat(str,"0");
        }

        snprintf(num, 4, "%d", i);

        strcat(str,num);
        strcat(str,".ppm");

        // On fait l'algorithme
        
        R=LoadPPM_rgb8matrix(str,&nrl,&nrh,&ncl,&nch);
        
        matrice_fond = bmatrix(nrl,nrh,ncl,nch);
        matrice_fond=LoadPGM_bmatrix(pathFirstImage,&nrl,&nrh,&ncl,&nch);

        matrice_temp = bmatrix(nrl,nrh,ncl,nch);
        matrice_difference = bmatrix(nrl,nrh,ncl,nch);

        matrice_temp = RGB_B(R,nrl,nrh,ncl,nch);
        
        //On fait la soutraction
        for(j = nrl;j < nrh;j++){
		    for(k = ncl;k < nch;k++){
                matrice_difference[j][k] = (byte)(abs((int)matrice_temp[j][k] - (int)matrice_fond[j][k]));
            }
        }
        
        //On traite l'image
        matrice_difference = binarisation(matrice_difference,nrl,nrh,ncl,nch,SEUIL);
        matrice_difference = ouverture(matrice_difference,nrl,nrh,ncl,nch,2);
        matrice_difference = fermeture(matrice_difference,nrl,nrh,ncl,nch,5);
   

        //On envoie dans un dossier
        memset (dest, 0, sizeof (dest));
        strcat(dest,pathfinal);
        strcat(dest,"image_difference");
        strcat(dest,num);
        strcat(dest,".pgm");    
        SavePGM_bmatrix(matrice_difference,nrl,nrh,ncl,nch,dest);
           
    }
}

/**
* Fonction permettant de définir une image de fond - b (Filtre median)
*/
void Exercice5_b(char *path,char *imagename,char *pathfinal,char *pathFirstImage){
    byte **matrice_temp;
    byte **matrice_difference;
    rgb8 **I;
    long nrh,nrl,nch,ncl;

    I=LoadPPM_rgb8matrix(pathFirstImage,&nrl,&nrh,&ncl,&nch);
    
    matrice_difference = bmatrix(nrl,nrh,ncl,nch);
    matrice_difference = RGB_B(I,nrl,nrh,ncl,nch);    

	int **tab = NULL;
	tab = InitTab(nrh*nch,256);

    rgb8 **R;
    int i,j,k;
    char str[50] = "\0";
    char num[4];
    char dest[50] = "\0";

	// On initialise le tableau à 0
	for(i = 0;i<nrh*nch;i++){
		for(j = 0; j < 256;j++){
			tab[i][j] = 0;
			//printf("i = %d && j = %d \n",i,j);
		}
	}
	printf("INITIALISATION DONE \n");

    for(i = 1;i < NB_IMAGE; i++){
        // On obtient le bon nom du fichier pour i :
        memset (str, 0, sizeof (str));
        strcat(str,path);
        strcat(str,imagename);
        if(i < 10){
            strcat(str,"00");
        }
        else if((i < 100)&&(i>=10)){
            strcat(str,"0");
        }

        snprintf(num, 4, "%d", i);

        strcat(str,num);
        strcat(str,".ppm");

        // On fait l'algorithme
        
        R=LoadPPM_rgb8matrix(str,&nrl,&nrh,&ncl,&nch);
        matrice_temp = bmatrix(nrl,nrh,ncl,nch);
        matrice_temp = RGB_B(R,nrl,nrh,ncl,nch);

        //printf("Start Moyenne glissante \n");
        //On fait la moyenne glissante
        for(j = nrl;j < nrh;j++){
		    for(k = ncl;k < nch;k++){
				tab[j*k][matrice_temp[j][k]]++;
            }
        }
  		//printf("etape %d \n",i);

    }
    printf("Moyenne glissante done\n");
    printf("Mediane");
    int flag = 0;
	int compteur = 0;
	int iter_ligne = 0;
	int iter_colonne = 0;
	//On fait ensuite la médiane
	for(i = 0;i<nrh*nch;i++){
        compteur = 0;
		iter_colonne++;
		if(iter_colonne == nch){
			iter_colonne = 0;
			iter_ligne++;
		} 
		for(j = 0; j < 256;j++){
			compteur = compteur + tab[i][j];
			if((compteur < NB_IMAGE/2)&&(flag == 0)){
				matrice_difference[iter_ligne][iter_colonne] = j;
                printf("i = %d && valeur = %d\n ",i,j);
                flag = 1;
			}
		}
	}
    printf("Mediane done");

    //On envoie dans un dossier
    memset (dest, 0, sizeof (dest));
    strcat(dest,pathfinal);
    strcat(dest,"mediane");
    strcat(dest,num);
    strcat(dest,".pgm");    
    SavePGM_bmatrix(matrice_difference,nrl,nrh,ncl,nch,dest);
    
}


/**
* Fonction permettant de définir une image de fond - a (Moyenne glissante)
*/
void Exercice5_a(char *path,char *imagename,char *pathfinal,char *pathFirstImage){
    byte **matrice_temp;
    byte **matrice_difference;
    rgb8 **I;
    long nrh,nrl,nch,ncl;
 
    I=LoadPPM_rgb8matrix(pathFirstImage,&nrl,&nrh,&ncl,&nch);
    
    matrice_difference = bmatrix(nrl,nrh,ncl,nch);
    matrice_difference = RGB_B(I,nrl,nrh,ncl,nch);    

    rgb8 **R;
    int i,j,k;
    char str[50] = "\0";
    char num[4];
    char dest[50] = "\0";
    for(i = 1;i < NB_IMAGE; i++){
        // On obtient le bon nom du fichier pour i :
        memset (str, 0, sizeof (str));
        strcat(str,path);
        strcat(str,imagename);
        if(i < 10){
            strcat(str,"00");
        }
        else if((i < 100)&&(i>=10)){
            strcat(str,"0");
        }

        snprintf(num, 4, "%d", i);

        strcat(str,num);
        strcat(str,".ppm");

        // On fait l'algorithme
        
        R=LoadPPM_rgb8matrix(str,&nrl,&nrh,&ncl,&nch);

        matrice_temp = bmatrix(nrl,nrh,ncl,nch);

        matrice_temp = RGB_B(R,nrl,nrh,ncl,nch);

        //On fait la moyenne glissante
        for(j = nrl;j < nrh;j++){
		    for(k = ncl;k < nch;k++){
                matrice_difference[j][k] = (matrice_temp[j][k] + matrice_difference[j][k])/2;
            }
        }
        
    }

    //On envoie dans un dossier
    memset (dest, 0, sizeof (dest));
    strcat(dest,pathfinal);
    strcat(dest,"moyenne_glissante");
    strcat(dest,num);
    strcat(dest,".pgm");    
    SavePGM_bmatrix(matrice_difference,nrl,nrh,ncl,nch,dest);
    
}

/**
* Fonction permettant de détecter le mouvement par différence d’images consécutives en fonction d'un nom
*/
void Exercice4(char *path,char *imagename,char *pathfinal){
    byte **matrice_temp;
    byte **matrice_temp_next;
    byte **matrice_difference;

    rgb8 **R;
    rgb8 **R_next;
    
    int i,j,k;
    int i_next;
    long nrh,nrl,nch,ncl;
    char str[50] = "\0";
    char str_next[50] = "\0";
    char num[4];
    char dest[50] = "\0";
    for(i = 1;i < 965; i++){
        // On obtient le bon nom du fichier pour i :
        memset (str, 0, sizeof (str));
        strcat(str,path);
        strcat(str,imagename);
        if(i < 10){
            strcat(str,"00");
        }
        else if((i < 100)&&(i>=10)){
            strcat(str,"0");
        }

        snprintf(num, 4, "%d", i);

        strcat(str,num);
        strcat(str,".ppm");

        // On obtient le bon nom du fichier pour i_next :
        i_next = i + 1 ; 
        memset (str_next, 0, sizeof (str_next));
        strcat(str_next,path);
        strcat(str_next,imagename);
        if(i_next < 10){
            strcat(str_next,"00");
        }
        else if((i_next < 100)&&(i_next>=10)){
            strcat(str_next,"0");
        }

        snprintf(num, 4, "%d", i_next);

        strcat(str_next,num);
        strcat(str_next,".ppm");

        // On fait l'algorithme
        
        R=LoadPPM_rgb8matrix(str,&nrl,&nrh,&ncl,&nch);
        R_next=LoadPPM_rgb8matrix(str_next,&nrl,&nrh,&ncl,&nch);

        matrice_temp = bmatrix(nrl,nrh,ncl,nch);
        matrice_temp_next = bmatrix(nrl,nrh,ncl,nch);
        matrice_difference = bmatrix(nrl,nrh,ncl,nch);

        matrice_temp = RGB_B(R,nrl,nrh,ncl,nch);
        matrice_temp_next = RGB_B(R_next,nrl,nrh,ncl,nch);
        
        //On fait la soutraction
        for(j = nrl;j < nrh;j++){
		    for(k = ncl;k < nch;k++){
                matrice_difference[j][k] = (byte)(abs((int)matrice_temp[j][k] - (int)matrice_temp_next[j][k]));  
            }
        }
        
        //On envoie dans un dossier
        memset (dest, 0, sizeof (dest));
        strcat(dest,pathfinal);
        strcat(dest,"image_difference");
        strcat(dest,num);
        strcat(dest,".ppm");    
        SavePGM_bmatrix(matrice_difference,nrl,nrh,ncl,nch,dest);
        printf("path %s\n = ",str_next);
           
    }
}

/**
* Fonction permettant de faire la fermeture
*/
byte **fermeture(byte **m,long nrl,long nrh,long ncl,long nch,int iteration){
	int i;
	for(i = 0;i < iteration;i++){
		m = dilatation(m,nrl,nrh,ncl,nch);
	}
	for(i = 0;i < iteration;i++){
		m = erosion(m,nrl,nrh,ncl,nch);
	}
	return m;
}

/**
* Fonction permettant de faire l'ouverture
*/
byte **ouverture(byte **m,long nrl,long nrh,long ncl,long nch,int iteration){
	int i;
	for(i = 0;i < iteration;i++){
		m = erosion(m,nrl,nrh,ncl,nch);
	}
	for(i = 0;i < iteration;i++){
		m = dilatation(m,nrl,nrh,ncl,nch);
	}
	return m;
}

/**
* Fonction permettant de faire l'érosion à l'aide 
*/
byte **erosion(byte **m,long nrl,long nrh,long ncl,long nch){
	long i,j;
	byte **matrice_temp;
	matrice_temp = bmatrix(nrl,nrh,ncl,nch);
	for(i = nrl + 1;i < nrh;i++){
		for(j = ncl + 1;j < nch;j++){
			if(m[i][j] == 255){
				if((m[i-1][j] == 0)||(m[i+1][j] == 0)||(m[i][j-1] == 0)||(m[i][j+1] == 0)||(m[i+1][j+1] == 0)||(m[i-1][j-1] == 0)||(m[i-1][j+1] == 0)||(m[i+1][j-1] == 0)){
					matrice_temp[i][j] = 0;
				}
				else{
					matrice_temp[i][j] = 255;
				}
			}
			else{
				matrice_temp[i][j] = 0;
			}
		}
	}
	return matrice_temp;
}


/**
* Fonction permettant de faire la dilation
*/
byte **dilatation(byte **m,long nrl,long nrh,long ncl,long nch){
	long i,j;
	byte **matrice_temp;
	matrice_temp = bmatrix(nrl,nrh,ncl,nch);
	for(i = nrl + 1;i < nrh;i++){
		for(j = ncl + 1;j < nch;j++){
			if(m[i][j] == 0){
				if((m[i-1][j] == 255)||(m[i+1][j] == 255)||(m[i][j-1] == 255)||(m[i][j+1] == 255)||(m[i+1][j+1] == 255)||(m[i-1][j+1] == 255)||(m[i-1][j-1] == 255)||(m[i+1][j-1] == 255)){
					matrice_temp[i][j] = 255;
				}
			}
			else{
				matrice_temp[i][j] = 255;
			}
		}
	}
	return matrice_temp;
}


/**
* Binarise l'image à partir d'un certain seuil passé en paramètre
*/
byte **binarisation(byte **m,long nrl,long nrh,long ncl,long nch,int seuil){
	long i,j;
	byte **matrice_temp;
	matrice_temp = bmatrix(nrl,nrh,ncl,nch);

	for(i = nrl;i < nrh;i++){
		for(j = ncl;j < nch;j++){
			if(m[i][j] < seuil)
				matrice_temp[i][j] = (byte)0;
			else
				matrice_temp[i][j] = (byte)255;
		}
	}
	return matrice_temp;
}

/**
* RGB ----> BINARY
*/
byte **RGB_B(rgb8 **m,long nrl,long nrh,long ncl,long nch){
	long i,j;
	byte **matrice_temp;
	matrice_temp = bmatrix(nrl,nrh,ncl,nch);
    int color;

	for(i = nrl;i < nrh;i++){
		for(j = ncl;j < nch;j++){
            color = (m[i][j].r + m[i][j].g + m[i][j].b)/3;
            matrice_temp[i][j] = color;
		}
	}
	return matrice_temp;
}

int main(void){

	
	long nrh,nrl,nch,ncl;
	rgb8 **I;
	byte **R;

	int i,j;
	//I=LoadPPM_rgb8matrix("img/lbox001.ppm",&nrl,&nrh,&ncl,&nch);

	//R = RGB_B(I,nrl,nrh,ncl,nch);

	//R = ouverture(R,nrl,nrh,ncl,nch,10);
    
    //Exercice6("Sequences/Lbox/ppm/\0","lbox\0","Sequences/Lbox/imgdone-mediane/","Sequences/Lbox/mediane869.pgm");
    Harris("Sequences/Fomd/ppm/\0","fomd\0","Sequences/Fomd/imgdone-harris/","Sequences/Fomd/ppm/fomd001.ppm",10);
    //Exercice5_b("Sequences/Fomd/ppm/\0","fomd\0","Sequences/Fomd/","Sequences/Fomd/ppm/fomd001.ppm");
    
	//SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"imgdone/test.pgm");

	// Creation d'une copie de l'image avec un nom different
	//SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"imgdone/carre_bina.pgm");
	//SavePGM_bmatrix(matrice_temp,nrl,nrh,ncl,nch,"imgdone/carre_bina.pgm");
	//free_bmatrix(R,nrl,nrh,ncl,nch);
	return 0;

}
