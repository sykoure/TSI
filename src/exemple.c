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

int SEUIL = 40;

byte **dilatation(byte **m,long nrl,long nrh,long ncl,long nch);
byte **erosion(byte **m,long nrl,long nrh,long ncl,long nch);
byte **RGB_B(rgb8 **m,long nrl,long nrh,long ncl,long nch);

/**byte passeBas(rgb8 **m,int i,int j){
	int ligne,colonne,compteur_l,compteur_c;
	int tempx = 0;
	int tempy = 0;
	int magnitude;
	compteur_c = 0;
	compteur_l = 0;
	for(ligne = i-1;ligne<=i+1;ligne++){
		for(colonne = j-1;colonne<=j+1;colonne++){
			tempx = tempx + xGradient[compteur_c][compteur_l] * m[ligne][colonne].r;
			tempy = tempy + yGradient[compteur_c][compteur_l] * m[ligne][colonne].r;
			if(compteur_c == 2){
				compteur_c = 0;
				compteur_l++;
			}
			else{
				
				compteur_c++;
			}
		}
	}
	if(tempx < 0){
		tempx = 0;
	}
	if(tempy < 0){
		tempy = 0;
	}

	magnitude = sqrt(tempx*tempx + tempy*tempy);	
	if(magnitude < SEUIL){
		magnitude = 0;
	}
	else{
		magnitude = 255;
	}
	printf("temp = %d && matrice = %d ; i = %d && j = %d\n",magnitude,m[i][j],i,j);
	return (byte)(magnitude);
}		


void convolution(rgb8 **m,long nrl,long nrh,long ncl,long nch){
	byte **matrice_convo;
	matrice_convo = bmatrix(nrl,nrh,ncl,nch);

	long i,j;
	for(i = nrl;i < nrh;i++){
		for(j = ncl;j < nch;j++){
			if((i == 0)||(j == 0)||(i == nrh-1)||(j == nch - 1)){
				matrice_convo[i][j] = 128;

			}
			else{
				matrice_convo[i][j] = passeBas(m,i,j);
			}
		}
	}
	SavePGM_bmatrix(matrice_convo,nrl,nrh,ncl,nch,"Seuil.pgm");
	free_bmatrix(matrice_convo,nrl,nrh,ncl,nch);
}

*/


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
        //printf("path %s\n = ",str_next);
           
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
    
    Exercice4("Sequences/Fomd/ppm/\0","fomd\0","Sequences/Fomd/imgdone/");
    
	//SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"imgdone/test.pgm");

	// Creation d'une copie de l'image avec un nom different
	//SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"imgdone/carre_bina.pgm");
	//SavePGM_bmatrix(matrice_temp,nrl,nrh,ncl,nch,"imgdone/carre_bina.pgm");
	//free_bmatrix(R,nrl,nrh,ncl,nch);

	

	return 0;

}
