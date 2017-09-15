//[[Rcpp::depends(RcppGSL)]]


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <algorithm>
#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_math.h> 
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <sstream> // pour compteColonnes
#include <math.h> // necessaire pour que floor soit defini 
#include <fstream> // lecture et ecriture de fichiers
#include <stdio.h>
#include <vector>

#include <Rcpp.h>
#include <RcppGSL.h>

using namespace std;
using namespace Rcpp;

//provient de MatriceEcarts.h



class MatriceEcarts
{
 public :

  // Constructeur : utilise les emplacements des points du planXP sous forme de gsl_matrix*.
  MatriceEcarts(gsl_matrix* const coord_points);

  // Nombre de points du planXP.
  unsigned int getNbPoints() const;
  // Dimension de l'espace dans lequel vivent les points du planXP.
  unsigned int getDimension() const;

  // 2 manieres differentes de fournir un ecart entre 2 points identifies par leur numero (ce qui permet de retrouver la bonne ligne dans *m_ecarts si numero_pt1!=numero_pt2 et de regarder *m_ecart_nul si numero_pt1==numero_pt2) :
  // 1) getEcartCopie cree un nouveau gsl_vector dans lequel il stocke ledit ecart ; --> possible amelioration : getEcartCopie pourrait prendre en parametre un gsl_vector* pour y inscrire les valeurs, et n'aurait donc pas a allouer des cases memoires a un nouveau gsl_vector.
  // 2) getEcartView renvoie un gsl_vector_view, et &RESULTAT.vector est un pointeur vers la ligne de *m_ecarts concernee. INUTILE DORENAVANT. DEVRAIT VRAISEMBLABLEMENT ETRE SUPPRIME.
  gsl_vector* getEcartCopie(unsigned int const& numero_pt1, unsigned int const& numero_pt2);
  //gsl_vector_view getEcartView(unsigned int const& numero_pt1, unsigned int const& numero_pt2);
  void getEcartCopie(gsl_vector* destination, unsigned int const& numero_pt1, unsigned int const& numero_pt2); // Dans cette version de getEcartCopie, l'ecart est copie dans le gsl_vector pris en parametre destination, et non renvoye. Attention, aucun test n'est fait pour savoir si destination a la taille adequate.

  ~MatriceEcarts();
   
 protected :

 unsigned int m_nb_points;
 unsigned int m_dimension; 

  gsl_matrix *m_ecarts;
  //gsl_vector *m_ecart_nul; // a priori INUTILE si getEcartView n'est pas utilise

};


//fin de MatriceEcarts.h


//provient de MatriceEcarts.cpp

//#include "MatriceEcarts.h"

MatriceEcarts::MatriceEcarts( gsl_matrix* const coord_points) : m_nb_points(coord_points->size1), m_dimension(coord_points->size2)
{
  // Cardinal de l'ensemble des ecarts differents de (0,0,...,0) entre points. 
  //Correspond au nombre de coefficients non diagonaux superieurs (ou inferieurs) d'une matrice carree d'ordre m_nb_points. 
  //Aussi egal au coefficient binomial d'indices m_nb_points et 2.
  unsigned int nb_lignes = m_nb_points * (m_nb_points-1) / 2; 
 
   
  m_ecarts = gsl_matrix_alloc(nb_lignes , m_dimension );


  // m_ecart_nul = gsl_vector_calloc(m_dimension); // allocation + initialisation a zero // INUTILE si getEcartView n'est pas utilise


  unsigned int numero_ligne(0);

  gsl_vector *point1 = gsl_vector_alloc(m_dimension);
  gsl_vector *point2 = gsl_vector_alloc(m_dimension);
  for (int i=0; i<m_nb_points; i++)
    {
      gsl_matrix_get_row(point1,coord_points, i) ; // Ligne i de coord_points copiee dans point1.
      for(int j=i+1; j<m_nb_points; j++)
  {
    gsl_matrix_get_row(point2, coord_points, j); // Ligne j de coord_points copiee dans point2
    gsl_vector_sub(point2, point1 ); // point2 = point2 - point1
    gsl_matrix_set_row( m_ecarts, numero_ligne, point2 ); // point2 (ex point2 - point1) copie dans la ligne numero_ligne de m_ecarts.

    // Noter que point1 n'a pas ete modifie au sein de cette 2e boucle for.
    numero_ligne++;
  }
    }

  // Liberation de la memoire.
  gsl_vector_free(point1);
  gsl_vector_free(point2);
}

MatriceEcarts::~MatriceEcarts()
{
  gsl_matrix_free(m_ecarts);
  // gsl_vector_free(m_ecart_nul); //m_ecart_nul INUTILE si getEcartView n'est pas utilise
}

unsigned int MatriceEcarts::getNbPoints() const
{
  return m_nb_points;
}

unsigned int MatriceEcarts::getDimension() const
{
  return m_dimension;
}

gsl_vector* MatriceEcarts::getEcartCopie(unsigned int const& numero_pt1, unsigned int const& numero_pt2)
{
  unsigned int numero_min ( min(numero_pt1,numero_pt2) );
  unsigned int numero_max ( max(numero_pt1,numero_pt2) );
  
  if(numero_max>=m_nb_points)
    {
      cout << "Erreur dans MatriceEcarts::getEcartCopie :" << endl;
      cout << "Il y a " << m_nb_points << " points, or vous voulez regarder un inexistant point numero " << numero_max << endl;
      cout << "NB : La numerotation commence a 0." << endl;
      cout << "Je renvoie le pointeur NULL." << endl << endl;
      return NULL;
    }
  
  else if(numero_min==numero_max)
    {
      gsl_vector *vecteur_ecart = gsl_vector_calloc(m_dimension); // vecteur de taille m_dimension rempli de zeros
      return vecteur_ecart;
    }

  else
    {
      gsl_vector *vecteur_ecart = gsl_vector_alloc(m_dimension);

      /* 
Les points de planXP sont numerotes (numero de ligne dans la gsl_matrix qui les contenait initialement -->cf coord_points dans le constructeur).
L'ecart entre 2 points peut donc etre indexe par 2 nombres correspondant aux numeros des 2 points. Ainsi, 0 11 indique l'ecart entre le point 0 (1er) et le point 11 (12e). Mais 11 0 represente le meme ecart (on rappelle que si l'ecart est multidimensionnel, il n'a pas de signe, meme si on n'a pas fait l'effort de les retirer dans le constructeur). Par convention, on designe toujours l'ecart par (plus petit numero) (plus grand numero) : 0 11 dans notre exemple.
A partir de la, on peut definir sur la base de la relation d'ordre entre numeros de points (le point 4 passe avant le point 11) une relation d'ordre entre ecarts, en definissant le successeur de a b par a b+1 s'il existe un point b+1, et a+1 a+2 sinon.
Les ecarts ont ete ranges dans cet ordre par le constructeur.
Du coup, on accede a l'ecart (numero_min) (numero_max) par :
 */

      gsl_matrix_get_row(vecteur_ecart, m_ecarts, numero_min*(m_nb_points-1) - numero_min*(numero_min-1)/2 + numero_max - 1 -  numero_min ); // Le -1 est du au fait que l'ecart entre points identique (0 0 ... 0) n'a ete stocke nulle part, alors qu'il serait logiquement le premier s'il l'avait ete.

      return vecteur_ecart;
    }
}



void MatriceEcarts::getEcartCopie(gsl_vector* destination, unsigned int const& numero_pt1, unsigned int const& numero_pt2)
{
  unsigned int numero_min ( min(numero_pt1,numero_pt2) );
  unsigned int numero_max ( max(numero_pt1,numero_pt2) );
  
  if(numero_max>=m_nb_points)
    {
      cout << "Erreur dans MatriceEcarts::getEcartCopie :" << endl;
      cout << "Il y a " << m_nb_points << " points, or vous voulez regarder un inexistant point numero " << numero_max << endl;
      cout << "NB : La numerotation commence a 0." << endl;
      cout << "Je ne fais rien." << endl << endl;
    }
  
  else if(numero_min==numero_max)
    {
      gsl_vector_set_zero(destination);
    }

  else
    {

      /* 
Les points de planXP sont numerotes (numero de ligne dans la gsl_matrix qui les contenait initialement -->cf coord_points dans le constructeur).
L'ecart entre 2 points peut donc etre indexe par 2 nombres correspondant aux numeros des 2 points. Ainsi, 0 11 indique l'ecart entre le point 0 (1er) et le point 11 (12e). Mais 11 0 represente le meme ecart (on rappelle que si l'ecart est multidimensionnel, il n'a pas de signe, meme si on n'a pas fait l'effort de les retirer dans le constructeur). Par convention, on designe toujours l'ecart par (plus petit numero) (plus grand numero) : 0 11 dans notre exemple.
A partir de la, on peut definir sur la base de la relation d'ordre entre numeros de points (le point 4 passe avant le point 11) une relation d'ordre entre ecarts, en definissant le successeur de a b par a b+1 s'il existe un point b+1, et a+1 a+2 sinon.
Les ecarts ont ete ranges dans cet ordre par le constructeur.
Du coup, on accede a l'ecart (numero_min) (numero_max) par :
 */

      gsl_matrix_get_row(destination, m_ecarts, numero_min*(m_nb_points-1) - numero_min*(numero_min-1)/2 + numero_max - 1 -  numero_min ); // Le -1 est du au fait que l'ecart entre points identique (0 0 ... 0) n'a ete stocke nulle part, alors qu'il serait logiquement le premier s'il l'avait ete.
    }
}

/*
gsl_vector_view MatriceEcarts::getEcartView(unsigned int const& numero_pt1, unsigned int const& numero_pt2)
{
  unsigned int numero_min ( min(numero_pt1,numero_pt2) );
  unsigned int numero_max ( max(numero_pt1,numero_pt2) );
  
  if(numero_max>=m_nb_points)
    {
      cout << "Erreur dans MatriceEcarts::getEcartView :" << endl;
      cout << "Il y a " << m_nb_points << " points, or vous voulez regarder un inexistant point numero " << numero_max << endl;
      cout << "NB : La numerotation commence a 0." << endl;
      //      cout << "Je renvoie le gsl_vector a une entree egale a 0." << endl << endl;
      //      cout << "Je ne renvoie rien." << endl << endl;
      cout << "Je ne sais pas quoi renvoyer." << endl << endl;
      //return gsl_vector_calloc(1);
    }
  
  else if(numero_min==numero_max)
    {
      return gsl_vector_subvector(m_ecart_nul,0,m_dimension);
    }

  else
    {
      return gsl_matrix_row(m_ecarts,  numero_min*(m_nb_points-1) - numero_min*(numero_min-1)/2 + numero_max - 1 -  numero_min );
    }
}
*/


//fin de MatriceEcarts.cpp

//provient de Matern.h




class NoyauMatern
{
 public :

  //NoyauMatern(double const& regularite); // Constructeur devenu inutile ; la valeur initiale de longueurs_correlation est fixee dans le programme et non renseignee en lignes de commande par l'utilisateur.
  NoyauMatern(double const& regularite, gsl_vector* longueurs_correlation);
  NoyauMatern(NoyauMatern const &noyau_a_copier);

  ~NoyauMatern();

  double getRegularite() const;
  int getDimension() const; // dimension de l'espace dans lequel vivent les points de planXP
  double getUneLongueurCorrelation(unsigned int const& indice_longueur_correlation_a_voir) const;
  gsl_vector* getLongueursCorrelation() const;
  void modifieUneLongueurCorrelation(unsigned int const& indice_longueur_correlation_a_modifier, double const& nouvelle_valeur);

  //void setLongueursCorrelation(gsl_vector *nouvelles_longueurs_correlation); // Cette methode servait a faire des tests. Normalement devenue inutile, mais gardee dans les commentaires par acquis de conscience. Il va de soi que cette methode serait DANGEREUSE dans un usage normal.

  double calCorrelation1d(double const& ecart_1d) const; // calcule la correlation obtenue avec un noyau de Matern entre 2 points de l'axe reel, et une longeur de correlation de 1.0.

  // CE QUI SUIT APPARTENAIT ANCIENNEMENT AUX CLASSES FILLES

  virtual double calCorrelation(gsl_vector* ecart_rd) const =0;
  virtual double calCorrelationDerivee(gsl_vector* ecart_rd, unsigned int const& indice_derive) =0;

 protected :

  double m_regularite;
  gsl_vector *m_longueurs_correlation; // tableau des longueurs de correlation
};


gsl_matrix* calMatriceCorrelation(MatriceEcarts* const matrice_ecarts, NoyauMatern* noyau);
gsl_matrix* calMatriceCorrelationDerivee(MatriceEcarts* const matrice_ecarts, unsigned int const& indice_derive, NoyauMatern* noyau);


double calDensitePrior(gsl_matrix* matrice_correlation, gsl_matrix* matrice_correlation_derivee, int const& nombre_observations_effectives); // bool const& berger, const gsl_matrix* tendance, gsl_vector* observations, const gsl_vector* observations_sauvegarde);
double calDensitePosterior(gsl_vector* observations, gsl_matrix* matrice_correlation_cholesky, double const& densite_prior);


// densitePosterior aggrege les 4 fonctions ci-dessus pour successivement calculer la matrice de correlation, sa derivee, puis les densites a priori et a posteriori.
double densitePosterior(gsl_vector* observations, unsigned int const& indice_derive, MatriceEcarts* matrice_ecarts, NoyauMatern* noyau, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives);//bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);


void calInjectionREML(gsl_matrix* tendance, gsl_matrix* injection_orthogonal_tendance); //, bool const& berger);
//void calProjectionBerger(const gsl_matrix* tendance, gsl_matrix* matrice_correlation_derivee, gsl_vector* observations, const gsl_matrix* matrice_correlation_cholesky, const gsl_vector* observations_sauvegarde);



//fin de Matern.h

//provient de Matern.cpp



/*
NoyauMatern::NoyauMatern(double const& regularite) :  m_regularite(regularite), m_longueurs_correlation(0)
{
  if(regularite<=0) 
    {
      cout << "La regularite doit etre strictement positive !" << endl;
      //      exit();
    }

  cout << "Renseigner le gsl_vector des longueurs de correlation :" << endl;

  //  m_longueurs_correlation = new Tableau;

  //OLD:
  // m_longueurs_correlation = new double[dimension];

  // for(int i=0; i<dimension; i++)
  //   {
  //     cout << "Longueur de correlation numero " << i+1 << " : ";
  //     cin >> m_longueurs_correlation[i];
  //     cout << endl;
  //   }

}
*/


NoyauMatern::NoyauMatern(double const& regularite, gsl_vector* longueurs_correlation) : m_regularite(regularite)
{
  if(regularite<=0) 
    {
      cout << "La regularite doit etre strictement positive !" << endl;
      //      exit();
    }

  m_longueurs_correlation = gsl_vector_alloc(longueurs_correlation->size);
  gsl_vector_memcpy ( m_longueurs_correlation,  longueurs_correlation);
}

NoyauMatern::NoyauMatern(NoyauMatern const &noyau_a_copier) : m_regularite(noyau_a_copier.m_regularite)
{
  m_longueurs_correlation = gsl_vector_alloc(noyau_a_copier.m_longueurs_correlation->size);
  gsl_vector_memcpy(m_longueurs_correlation, noyau_a_copier.m_longueurs_correlation);
}


NoyauMatern::~NoyauMatern()
{
  gsl_vector_free(m_longueurs_correlation);
}

double NoyauMatern::getRegularite() const
{
  return m_regularite;
}

int NoyauMatern::getDimension() const
{
  return m_longueurs_correlation->size;
}

double NoyauMatern::getUneLongueurCorrelation(unsigned int const& indice_longueur_correlation_a_voir) const
{
  return gsl_vector_get(m_longueurs_correlation,indice_longueur_correlation_a_voir);
}

gsl_vector* NoyauMatern::getLongueursCorrelation() const
{
  gsl_vector* copie_longueurs_correlation = gsl_vector_alloc(m_longueurs_correlation->size);
  gsl_vector_memcpy(copie_longueurs_correlation, m_longueurs_correlation);
  return copie_longueurs_correlation;

  // // if(copie_dimension != m_dimension)
  // //   {
  // //     cout << "Probleme dans getLongueursCorrelation :" << endl;
  // //     cout << "Dimension donnee : " << copie_dimension << endl;
  // //     cout << "Dimension du NoyauMatern : " << m_dimension << endl;
  // //     cout << "Je sors immediatement de getLongueursCorrelation." << endl;
  // //   }

  // // else
  // //   {
  // //     for(int i=0; i<copie_dimension; i++)
  // //   {
  // //     copie_longueurs_correlation[i] = m_longueurs_correlation[i];
  // //   }
  // //   }

}

void NoyauMatern::modifieUneLongueurCorrelation(unsigned int const& indice_longueur_correlation_a_modifier, double const& nouvelle_valeur)
{
  gsl_vector_set(m_longueurs_correlation,indice_longueur_correlation_a_modifier,nouvelle_valeur);
}

/*
void NoyauMatern::setLongueursCorrelation(gsl_vector *nouvelles_longueurs_correlation)
{
  gsl_vector_memcpy(m_longueurs_correlation, nouvelles_longueurs_correlation);
  // for(int i=0; i< getDimension();i++)
  //   {
  //     m_longueurs_correlation[i] = nouvelles_longueurs_correlation[i];
  //   }
}
*/

double NoyauMatern::calCorrelation1d(double const& ecart_1d) const
{
  if(ecart_1d==0) // La correlation est 1 si l'ecart est nul.
    {
      return 1;
    }

  double num_standard( 2.0* sqrt(m_regularite) * fabs(ecart_1d) ); // num est "standard" car les longueurs de correlation sont supposees avoir deja ete prises en compte

  // cout << "num_standard = " << num_standard << endl;
  //cout<<"I"<<endl;

  if(num_standard>705.0) return 0.0; // Au-dela d'un certaine valeur, gsl_sf_bessel_Knu est en overflow. De toute maniere, la correlation tend vers 0 quand num_standard tend vers l'infini, et 700.0 est une tres grande valeur de num_standard

  //cout << gsl_sf_exp(m_regularite  * gsl_sf_log(num_standard))<<endl;
  // cout<<"II"<<endl;
  //cout<<"num_standard = "<<num_standard<<endl;
  //cout << gsl_sf_bessel_Knu(m_regularite,num_standard)<<endl;
  //cout<<"III"<<endl;
  //cout << " * = " <<  gsl_sf_exp(m_regularite  * gsl_sf_log(num_standard)) * gsl_sf_bessel_Knu(m_regularite,num_standard)<<endl;
  //cout<<"IV"<<endl;

      return gsl_sf_gammainv (m_regularite) / gsl_sf_exp((m_regularite -1)* M_LN2) * gsl_sf_exp(m_regularite  * gsl_sf_log(num_standard)) * gsl_sf_bessel_Knu(m_regularite,num_standard);

  //  return gsl_sf_gammainv (m_regularite) / gsl_sf_exp((m_regularite -1)* M_LN2) *pow(num_standard, m_regularite) * gsl_sf_bessel_Knu(m_regularite,num_standard);
}





gsl_matrix* calMatriceCorrelation(MatriceEcarts* const matrice_ecarts, NoyauMatern* noyau, const gsl_matrix* injection_orthogonal_tendance)
{
  // 2 matrices necessaires : matrice_correlation (retournee a la fin) et matrice_auxiliaire
  gsl_matrix* matrice_correlation = gsl_matrix_alloc(injection_orthogonal_tendance->size1,injection_orthogonal_tendance->size1);
  gsl_matrix* matrice_auxiliaire = gsl_matrix_alloc(matrice_ecarts->getNbPoints(), matrice_ecarts->getNbPoints());

  

  // matrice_correlation est d'abord remplie de zeros...
  gsl_matrix_set_zero(matrice_correlation);


  // cout<<"----------------"<<endl;
  // cout<<matrice_correlation->size1<<" "<<matrice_correlation->size2<<"   "<<matrice_auxiliaire->size1<<" "<<matrice_auxiliaire->size2<<endl;
  // cout<<"----------------"<<endl;

  // gsl_vector_view ecart_vue;

  //cout<<"Dans calMatriceCorrelation."<<endl;
  //cout<<"Les longueurs_correlation sont " << noyau->getUneLongueurCorrelation(0) << " et " << noyau->getUneLongueurCorrelation(1) << endl << endl << endl;

   
  // ... puis on calcule ses coefficients triangulaires superieurs.

  for(int i=0;i<matrice_correlation->size1;i++)
    {
      for(int j=i+1;j<matrice_correlation->size2;j++) // en principe size2==size1
  {
    gsl_vector* ecart_vue = matrice_ecarts->getEcartCopie(i,j); // Le nom ecart_vue est un peu mal choisi : il s'agit bien d'une copie de la ligne correspondante de la matrice_ecarts, et non d'un gsl_vector_view, qui donnerait acces a ladite ligne de matrice_ecarts.
    gsl_matrix_set(matrice_correlation,i,j,noyau->calCorrelation(ecart_vue));
    gsl_vector_free(ecart_vue);
  }
    }

  // matrice_auxiliaire devient la transposee de matrice_correlation
  gsl_matrix_transpose_memcpy(matrice_auxiliaire,matrice_correlation);

  // Les coefficients triangulaires inferieurs de matrice_correlation sont completes par addition de sa transposee matrice_auxiliaire.
  gsl_matrix_add(matrice_correlation,matrice_auxiliaire);

  // matrice_auxiliaire est redefinie comme la matrice identite.
  gsl_matrix_set_identity(matrice_auxiliaire);

  // L'addition de ces coefficients diagonaux (tous egaux a 1.0) acheve la creation de matrice_correlation
  gsl_matrix_add(matrice_correlation,matrice_auxiliaire);

  // liberation de la memoire
  gsl_matrix_free(matrice_auxiliaire);


         // FILE * f = fopen ("matrice_correlation.txt", "w");
         // gsl_matrix_fprintf (f, matrice_correlation,"%f");
         // fclose (f);


  
  //NOUVEAU : il faut restreindre matrice_correlation a l'orthogonal de tendance.

  //La restriction est plus petite : elle a meme dimension que l'espace engendre par l'orthogonal a tendance (nb de colonnes de injection_orthogonal_tendance).
  gsl_matrix* matrice_correlation_orthogonal_tendance = gsl_matrix_alloc(injection_orthogonal_tendance->size2, injection_orthogonal_tendance->size2);
  
  gsl_matrix* auxiliaire2 = gsl_matrix_calloc(matrice_correlation->size1, injection_orthogonal_tendance->size2); //matrice stockant un produit matriciel intermediaire

  gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,matrice_correlation,injection_orthogonal_tendance,0.0,auxiliaire2);// auxiliaire2 = matrice_correlation injection_orthogonal_tendance
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,injection_orthogonal_tendance,auxiliaire2,0.0,matrice_correlation_orthogonal_tendance);// matrice_correlation_orthogonal_tendance = injection_orthogonal_tendance^T auxiliaire2

  gsl_matrix_free(auxiliaire2);
  gsl_matrix_free(matrice_correlation);

  //gsl_matrix* echange = matrice_correlation;// echange est un pointeur vers matrice_correlation, pour eviter une fuite de memoire
  //matrice_correlation = matrice_correlation_orthogonal_tendance;// matrice_correlation pointe desormais vers matrice_correlation_orthogonal_tendance
  //gsl_matrix_free(echange);// liberation de la memoire de echange (ex matrice_correlation)


  // {//cout << matrice_correlation_orthogonal_tendance->size1 << matrice_correlation_orthogonal_tendance->size2 << endl;  
  //       FILE * g = fopen ("matrice_correlation_orthogonal_tendance.txt", "w");
  //       gsl_matrix_fprintf (g, matrice_correlation_orthogonal_tendance,"%f");
  //       fclose (g);
  //     }
  //cout << "L'adresse de matrice_correlation_orthogonal_tendance est : " << matrice_correlation_orthogonal_tendance << endl;


    return matrice_correlation_orthogonal_tendance;

}

gsl_matrix* calMatriceCorrelationDerivee(MatriceEcarts* const matrice_ecarts, unsigned int const& indice_derive, NoyauMatern* noyau, const gsl_matrix* injection_orthogonal_tendance)
{
  // 2 matrices necessaires : matrice_correlation_derivee (retournee a la fin) et matrice_auxiliaire
  gsl_matrix* matrice_correlation_derivee = gsl_matrix_alloc(injection_orthogonal_tendance->size1, injection_orthogonal_tendance->size1);
  gsl_matrix* matrice_auxiliaire = gsl_matrix_alloc(matrice_ecarts->getNbPoints(), matrice_ecarts->getNbPoints());

  // matrice_correlation_derivee est d'abord remplie de zeros...
  gsl_matrix_set_zero(matrice_correlation_derivee);

    //gsl_vector_view ecart_vue;
  //cout<<"Dans calMatriceCorrelationDerivee."<<endl;
  //cout<<"Les longueurs_correlation sont " << noyau->getUneLongueurCorrelation(0) << " et " << noyau->getUneLongueurCorrelation(1) << endl << endl << endl;
 
  // on calcule ses coefficients triangulaires superieurs.

  gsl_vector* ecart_vue = gsl_vector_alloc(matrice_ecarts->getDimension());
  
  for(int i=0;i<matrice_correlation_derivee->size1;i++)
    {
      for(int j=i+1;j<matrice_correlation_derivee->size2;j++) // en principe size2==size1
  {
    matrice_ecarts->getEcartCopie(ecart_vue,i,j);  // Le nom ecart_vue est un peu mal choisi : il s'agit bien d'une copie de la ligne correspondante de la matrice_ecarts, et non d'un gsl_vector_view, qui donnerait acces a ladite ligne de matrice_ecarts.
    
    //cout<<i<<" "<<j<<endl;
    gsl_matrix_set(matrice_correlation_derivee,i,j,noyau->calCorrelationDerivee(ecart_vue,indice_derive));

    //if(i==0) cout<<gsl_matrix_get(matrice_correlation_derivee,i,j)<<endl;
    //if(i==0) cout<<noyau->calCorrelationDerivee(matrice_ecarts->getEcartCopie(i,j),indice_derive) << endl;
    //if(i==0) cout<<"---"<<endl;
  }
    }
  
  gsl_vector_free(ecart_vue);

  //cout << "minimum de matrice_correlation_derivee avant fin de la construction = " <<gsl_matrix_min(matrice_correlation_derivee)<<" et maximum = " << gsl_matrix_max(matrice_correlation_derivee)<<endl;

  // matrice_auxiliaire devient la transposee de matrice_correlation_derivee
  gsl_matrix_transpose_memcpy(matrice_auxiliaire,matrice_correlation_derivee);

  // Les coefficients triangulaires inferieurs de matrice_correlation_derivee sont completes par addition a matrice_correlation_derivee de sa transposee matrice_auxiliaire.
  gsl_matrix_add(matrice_correlation_derivee,matrice_auxiliaire);

  // Pas besoin de s'occuper des coefficients diagonaux de matrice_correlation_derivee, puisqu'ils sont de toute maniere nuls


  //cout << "minimum de matrice_correlation_derivee = " <<gsl_matrix_min(matrice_correlation_derivee)<<" et maximum = " << gsl_matrix_max(matrice_correlation_derivee)<<endl;

  
  // Liberation de la memoire
  gsl_matrix_free(matrice_auxiliaire);


         // FILE * f = fopen ("matrice_correlation_derivee.txt", "w");
         // gsl_matrix_fprintf (f, matrice_correlation_derivee,"%f");
         // fclose (f);


  
  //NOUVEAU : il faut restreindre matrice_correlation_derivee a l'orthogonal de tendance.

  //La restriction est plus petite : elle a meme dimension que l'espace engendre par l'orthogonal a tendance (nb de colonnes de injection_orthogonal_tendance).
  gsl_matrix* matrice_correlation_derivee_orthogonal_tendance = gsl_matrix_alloc(injection_orthogonal_tendance->size2, injection_orthogonal_tendance->size2);
  
  gsl_matrix* auxiliaire2 = gsl_matrix_calloc(matrice_correlation_derivee->size1, injection_orthogonal_tendance->size2); //matrice stockant un produit matriciel intermediaire

  gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,matrice_correlation_derivee,injection_orthogonal_tendance,0.0,auxiliaire2);// auxiliaire2 = matrice_correlation_derivee injection_orthogonal_tendance
  
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,injection_orthogonal_tendance,auxiliaire2,0.0,matrice_correlation_derivee_orthogonal_tendance);// matrice_correlation_derivee_orthogonal_tendance = injection_orthogonal_tendance^T auxiliaire2

  gsl_matrix_free(auxiliaire2);
  gsl_matrix_free(matrice_correlation_derivee);

  //gsl_matrix* echange = matrice_correlation_derivee;// echange est un pointeur vers matrice_correlation, pour eviter une fuite de memoire
  // matrice_correlation_derivee = matrice_correlation_derivee_orthogonal_tendance;// matrice_correlation pointe desormais vers matrice_correlation_orthogonal_tendance
  //gsl_matrix_free(echange);// liberation de la memoire de echange (ex matrice_correlation)


         // FILE * g = fopen ("matrice_correlation_derivee_orthogonal_tendance.txt", "w");
         // gsl_matrix_fprintf (g, matrice_correlation_derivee_orthogonal_tendance,"%f");
         // fclose (g);

  
  return matrice_correlation_derivee_orthogonal_tendance;
  
}










double calDensitePrior(gsl_matrix* matrice_correlation, gsl_matrix* matrice_correlation_derivee, int const& nombre_observations_effectives ) //bool const& berger, const gsl_matrix* tendance, gsl_vector* observations, const gsl_vector* observations_sauvegarde)
{


  // gsl_matrix* auxiliaire = gsl_matrix_calloc(matrice_correlation_derivee->size1, matrice_correlation_derivee->size2);  
  // gsl_blas_dsymm(CblasLeft, CblasLower, 1.0, matrice_correlation_derivee, projection, 0.0, auxiliaire); // auxiliaire = matrice_correlation_derivee %*% projection
  // gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,projection,auxiliaire,0.0,matrice_correlation_derivee); // matrice_correlation_derivee = projection %*% auxiliaire
  // gsl_matrix_free(auxiliaire);

  
  //gsl_matrix* matrice_correlation = calMatriceCorrelation(matrice_ecarts, noyau);

  // gsl_matrix* matrice_correlation_derivee = calMatriceCorrelationDerivee(matrice_ecarts, indice_derive, noyau);

  // // Decoposition LU de matrice_correlation, afin de pouvoir calculer son determinant et son inverse
  // gsl_permutation* perm = gsl_permutation_alloc(matrice_correlation->size1);
  // int* sig = new int(0);
  // gsl_linalg_LU_decomp(matrice_correlation,perm,sig);

  // double det = gsl_linalg_LU_det(matrice_correlation,sig);

  // gsl_linalg_LU_invert(matrice_correlation,perm,sig);

  // gsl_blas_dtrsm (CblasLeft, CblasLower, CblasNoTrans, CblasUnit,1.0, matrice_correlation, matrice_correlation_derivee);

  // gsl_blas_dtrsm (CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,1.0, matrice_correlation, matrice_correlation_derivee);

  // gsl_vector_view diag_simple = gsl_matrix_diagonal(matrice_correlation_derivee);

  // double trace_simple = gsl_blas_dasum(diag_simple);


  //Decomposition de Cholesky de matrice_correlation, afin de pouvoir calculer son determinant et son inverse.
  gsl_linalg_cholesky_decomp(matrice_correlation);

  // Projection liee au prior de Berger s'il y a lieu PLUS NECESSAIRE !!!
  //if(berger) calProjectionBerger(tendance,matrice_correlation_derivee,observations,matrice_correlation,observations_sauvegarde);

  // matrice_correlation_derivee= L^{-1} matrice_correlation_derivee, avec L la transformee de Cholesky de la matrice de correlation (ie L L^{T} est la matrice de correlation). 
  // Attention, du fait de l'operation ci-dessus, matrice_correlation est maintenant == L.
  gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, matrice_correlation, matrice_correlation_derivee);

  // matrice_correlation_derivee = matrice_correlation_derivee L^{-T}, avec L comme ci-dessus et matrice_correlation_derivee == L^{-1} (matrice_correlation_derivee ancienne).
  // De ce fait, matrice_correlation_derivee nouvelle == L^{-1} matrice_correlation_derivee ancienne L^{-T}.
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1.0, matrice_correlation, matrice_correlation_derivee);

  // contrib_simple = -1.0/nb de points * Trace(matrice_correlation_derivee nouvelle).
  // Rappel : Trace(matrice_correlation_derivee nouvelle) == Trace (L^{-1} matrice_correlation_derivee ancienne L^{-T}) == Trace(matrice_correlation_derivee nouvelle L^{-T} L^{-1}) == Trace (matrice_correlation_derivee ancienne matrice_correlation ancienne^{-1}).
  gsl_vector_view diag_simple = gsl_matrix_diagonal(matrice_correlation_derivee);
  double contrib_simple = -1.0/nombre_observations_effectives * gsl_pow_2(gsl_blas_dasum(&diag_simple.vector)) ;

  // carre est matrice_correlation_derivee nouvelle prise au carre.
   gsl_matrix* carre = gsl_matrix_calloc(matrice_correlation_derivee->size1, matrice_correlation_derivee->size2); // En principe, size1 == size2
   gsl_blas_dsymm(CblasLeft,CblasLower,1.0,matrice_correlation_derivee,matrice_correlation_derivee,0.0,carre); // carre = 1.0 * matrice_correlation_derivee nouvelle matrice_correlation_derivee nouvelle + 0.0 * carre
   // Attention, si on essaie de remplacer carre dans les arguments de la fonction ci-dessus par matrice_correlation_derivee, on aboutit a la gsl_matrix rempie de zeros.
   // Remarque : matrice_correlation_derivee apparait 2 fois dans les arguments et est a chaque fois declaree const. Evidemment, carre n'est pas declaree const.


  //cout<<"minimum du simple = "<<gsl_matrix_min(matrice_correlation_derivee)<<"  maximum du simple = "<<gsl_matrix_max(matrice_correlation_derivee)<<endl;
  //cout<<"minimum du carre = "<<gsl_matrix_min(carre)<<"  maximum du carre = "<<gsl_matrix_max(carre)<<endl;

   // contrib_carre = Trace(matrice_correlation_derivee nouvelle matrice_correlation_derivee nouvelle) 
  gsl_vector_view diag_carre = gsl_matrix_diagonal(carre);
  //gsl_vector_view diag_carre = gsl_matrix_diagonal(matrice_correlation_derivee);
  double contrib_carre = gsl_blas_dasum(&diag_carre.vector);

  // Liberation de la memoire.
  gsl_matrix_free(carre);

  double densite_prior = sqrt( contrib_carre + contrib_simple );
  //cout<<"contrib_carre = "<<contrib_carre<<endl;;
  //cout<<"contrib_simple = "<<contrib_simple<<endl;;
  //cout<<"densite_prior = "<<densite_prior<<endl;;

  return densite_prior;
}

double calDensitePosterior(gsl_vector* observations, gsl_matrix* matrice_correlation_cholesky, double const& densite_prior)
{
  // observations_decorrelees est cree et on lui affecte observations.
  gsl_vector* observations_decorrelees = gsl_vector_alloc(observations->size);
  gsl_vector_memcpy(observations_decorrelees,observations);

  //cout<<"norme2 de observations_decorrelees avant la decorrelation = "<<gsl_blas_dnrm2(observations_decorrelees)<<endl;

  // det_matrice_correlation_cholesky ==  sqrt(determinant(matrice de correlation -- la vraie)) == determinant(matrice_correlation_cholesky))
  double det_matrice_correlation_cholesky = gsl_linalg_LU_det(matrice_correlation_cholesky,1); // gsl_linalg_LU_det calcule juste le produit des coefficients diagonaux de la matrice qui lui est soumise, puis multiplie par 1 ou -1, selon la valeur de l'entier pointe par le 2e argument.


  //cout << "det_matrice_correlation_cholesky = "<<det_matrice_correlation_cholesky<<endl;

 // observations_decorrelees = matrice_correlation_cholesky^{-1} observations
  gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, matrice_correlation_cholesky, observations_decorrelees); 

  //cout<<"norme2 de observations_decorrelees apres la decorrelation = "<<gsl_blas_dnrm2(observations_decorrelees)<<endl;

  double norme2 = gsl_blas_dnrm2 (observations_decorrelees); // norme2 == sqrt(observations^T matrice_correlation^{-1} observations)

  //cout<<"norme2 = "<<norme2<<endl;
  //cout<<"norme2^{"<<observations->size<<"} = "<<  gsl_pow_int(norme2, observations->size)<<endl;

  // Liberation de la memoire
  gsl_vector_free(observations_decorrelees);

  double densite_posterior = densite_prior / det_matrice_correlation_cholesky / gsl_pow_int(norme2, observations->size);
  //cout<<"vraisemblance = "<< 1.0 / det_matrice_correlation_cholesky / gsl_pow_int(norme2, observations->size)<<endl;
  //cout<<"densite_posterior = "<<densite_posterior;

  return densite_posterior;
}





double densitePosterior(gsl_vector* observations, unsigned int const& indice_derive, MatriceEcarts* matrice_ecarts, NoyauMatern* noyau, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{
  //gsl_matrix* matrice_correlation = gsl_matrix_alloc(observations->size, observations->size);
  //gsl_matrix* matrice_correlation_derivee = gsl_matrix_alloc(observations->size, observations->size);

  //cout<<"Marqueur5"<<endl;
  gsl_matrix* matrice_correlation = calMatriceCorrelation(matrice_ecarts, noyau,injection_orthogonal_tendance);

  //   cout << "L'adresse de matrice_correlation est maintenant: " << matrice_correlation << endl;

  
  // {cout << "Dimension de matrice_correlation3 : " << matrice_correlation->size1 <<"  " << matrice_correlation->size2 << endl;  
  //     FILE * f = fopen ("matrice_correlation3.txt", "w");
  //     gsl_matrix_fprintf (f, matrice_correlation,"%f");
  //     fclose (f);
  //   }
  //cout<<"Marqueur6"<<endl;
  gsl_matrix* matrice_correlation_derivee = calMatriceCorrelationDerivee(matrice_ecarts, indice_derive, noyau,injection_orthogonal_tendance);
  //cout<<"Marqueur7"<<endl;

  //cout << "Dim de matrice_correlation_derivee :" << matrice_correlation_derivee->size1 << matrice_correlation_derivee->size2 << endl;



  
  double densite_prior = calDensitePrior(matrice_correlation, matrice_correlation_derivee,nombre_observations_effectives); // berger,tendance,observations,observations_sauvegarde);
  //cout<<"Marqueur8"<<endl;
  // ATTENTION ! calDensitePrior a eu pour "effets secondaires" de :
  // 1) transformer matrice_correlation en sa decomposition de Cholesky (matrice triangulaire inferieure L telle que matrice_correlation de depart == L L^T.
  // 2) transformer matrice_correlation_derivee en ( L^{-1} matrice_correlation_derivee L^{-T} )^2 , avec L la matrice decrite ci-dessus (qui est designee dans la suite de cette fonction par "matrice_correlation")
  double densite_posterior = calDensitePosterior(observations, matrice_correlation, densite_prior);
  //cout<<"Marqueur9"<<endl;
  gsl_matrix_free(matrice_correlation);
  gsl_matrix_free(matrice_correlation_derivee);

  return densite_posterior;
}




void calInjectionREML(gsl_matrix* tendance, gsl_matrix* injection_orthogonal_tendance) //, bool const& berger)
{
  
  gsl_matrix_set_identity(injection_orthogonal_tendance);

  if(tendance!=NULL) // & !berger)
    {
      gsl_matrix* matriceR = gsl_matrix_alloc(tendance->size1, tendance->size2);
  gsl_matrix* auxiliaire = gsl_matrix_alloc(tendance->size1, tendance->size1);
  //gsl_matrix* auxiliaire2 = gsl_matrix_calloc(tendance->size1, tendance->size1);
      
  //gsl_matrix* r = gsl_matrix_alloc(tendance->size1, tendance->size2);
      gsl_vector* tau = gsl_vector_alloc(tendance->size2);
      //gsl_permutation* p = gsl_permutation_alloc(tendance->size2);
      //int* signum = new int(0);
      // gsl_vector* norm = gsl_vector_alloc(tendance->size2);

  //auxiliaire devient la matrice dont les n-p premieres colonnes sont la matrice P de memes dimensions que tendance (H) mais dont les colonnes sont deux a deux orthogonales
  gsl_linalg_QR_decomp(tendance,tau);
  gsl_linalg_QR_unpack(tendance,tau,auxiliaire,matriceR);
    //gsl_matrix_free(r);
  gsl_vector_free(tau);
  //gsl_permutation_free(p);
  //gsl_vector_free(norm);

  gsl_matrix_free(matriceR);

  gsl_matrix_memcpy(injection_orthogonal_tendance, &gsl_matrix_submatrix(auxiliaire,0,tendance->size2, tendance->size1, tendance->size1 - tendance->size2).matrix);

  gsl_matrix_free(auxiliaire);


   // gsl_matrix_view injection = gsl_matrix_submatrix(auxiliaire,0,0,tendance->size1,tendance->size2);


   // FILE * h = fopen ("injection.txt", "w");
   // gsl_matrix_fprintf (h, &injection.matrix,"%f");
   // fclose (h);

   //    FILE * f = fopen ("injection_orthogonal_tendance.txt", "w");
   // gsl_matrix_fprintf (f, injection_orthogonal_tendance,"%f");
   // fclose (f);
  
 //  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,&injection.matrix,&injection.matrix,0.0,auxiliaire2);


 //  // {cout << auxiliaire2->size1 << auxiliaire2->size2 << endl;  
 //  //   FILE * g = fopen ("auxiliaire2.txt", "w");
 //  //      gsl_matrix_fprintf (g, auxiliaire2,"%f");
 //  //      fclose (g);
 //  //    }
  
  
 //  gsl_matrix_free(auxiliaire);


 //  gsl_matrix_scale(auxiliaire2,-1); //(-1)*auxiliaire
 //  gsl_matrix* projection_orthogonal_tendance = gsl_matrix_alloc(tendance->size1, tendance->size1);
 //  gsl_matrix_set_identity(projection_orthogonal_tendance);
 //  gsl_matrix_add(projection_orthogonal_tendance,auxiliaire2);//projection_orthogonal_tendance:=identite-auxiliaire
 //  //cout << gsl_matrix_max(projection_orthogonal_tendance)<<endl;
 //  gsl_matrix_free(auxiliaire2);

 //  //Retourne la matrix_view des tendance->size2 premieres colonnes de projection_orthogonal_tendance (ie la matrice W).
 //  //return gsl_matrix_submatrix(projection_orthogonal_tendance,0,0,tendance->size1, tendance->size1 - tendance->size2);


 //  // {cout << projection_orthogonal_tendance->size1 << projection_orthogonal_tendance->size2 << endl;  
 //  //      FILE * f = fopen ("projection_orthogonal_tendance.txt", "w");
 //  //      gsl_matrix_fprintf (f, projection_orthogonal_tendance,"%f");
 //  //      fclose (f);
 //  //    }

 // //les n-p premieres colonnes de auxiliaire3 forment la matrice injection_orthogonal_tendance
 //  gsl_matrix* auxiliaire3 = gsl_matrix_alloc(tendance->size1, tendance->size1);
 //  r = gsl_matrix_alloc(tendance->size1,tendance->size1);
 //  tau = gsl_vector_alloc(tendance->size1);
 //  p = gsl_permutation_alloc(tendance->size1);
 //  norm = gsl_vector_alloc(tendance->size1);
 //  gsl_linalg_QRPT_decomp2(projection_orthogonal_tendance,auxiliaire3,r,tau,p,signum,norm);
 //  gsl_matrix_free(r);
 //  gsl_vector_free(tau);
 //  gsl_permutation_free(p);
 //  gsl_vector_free(norm);
 //  delete signum;

 //  gsl_matrix_view auxiliaire3_restreint = gsl_matrix_submatrix(auxiliaire3,0,0,tendance->size1, tendance->size1 -  tendance->size2);

 //  gsl_matrix_memcpy(injection_orthogonal_tendance,&auxiliaire3_restreint.matrix);

 //  gsl_matrix_free(auxiliaire3);

 //        // FILE * h = fopen ("injection_orthogonal_tendance.txt", "w");
 //        // gsl_matrix_fprintf (h, injection_orthogonal_tendance,"%f");
 //        // fclose (h);
  

  cout<<"Fin de la creation de l'injection"<<endl;
    }
}


// void calProjectionBerger(const gsl_matrix* tendance, gsl_matrix* matrice_correlation_derivee, gsl_vector* observations, const gsl_matrix* matrice_correlation_cholesky, const gsl_vector* observations_sauvegarde)
// {
//   gsl_matrix* matrice_correlation_inverse = gsl_matrix_alloc(matrice_correlation_cholesky->size1, matrice_correlation_cholesky->size2); // matrice_correlation_cholesky = L^(-1) (RAPPEL)
//   gsl_matrix_set_identity(matrice_correlation_inverse);
//   gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,matrice_correlation_cholesky,matrice_correlation_inverse); // matrice_correlation_inverse = L^(-1)
//   gsl_blas_dtrsm(CblasLeft,CblasLower,CblasTrans,CblasNonUnit,1.0,matrice_correlation_cholesky,matrice_correlation_inverse); // matrice_correlation_cholesky = l^(-T) L^(-1)
  
//   gsl_matrix* matrice_tendance_correlation_inverse = gsl_matrix_calloc(tendance->size2, tendance->size1);
//   gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,tendance,matrice_correlation_inverse,0.0,matrice_tendance_correlation_inverse); // matrice_tendance_correlation_inverse = tendance^T matrice_correlation_inverse

//   gsl_matrix_set_zero(matrice_correlation_inverse); // on va reinitialiser cette matrice en la remplissant de zeros
//   gsl_matrix* matrice_correlation_derivee_copie = matrice_correlation_inverse; // nouveau pointeur vers la matrice venant d'etre reinitialisee 
//   gsl_matrix_memcpy(matrice_correlation_derivee_copie,matrice_correlation_derivee); // copie de matrice_correlation_derivee dans matrice_correlation_derivee_copie

//   gsl_matrix* matrice_tendance_correlation_inverse_tendance = gsl_matrix_calloc(tendance->size2, tendance->size2);
//   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,matrice_tendance_correlation_inverse, tendance,0.0,matrice_tendance_correlation_inverse_tendance); // matrice_tendance_correlation_inverse_tendance = matrice_tendance_correlation_inverse tendance
//   //Decomposition de Cholesky de matrice_tendance_correlation_inverse_tendance, afin de pouvoir calculer son determinant et son inverse.
//   gsl_linalg_cholesky_decomp(matrice_tendance_correlation_inverse_tendance); //matrice_tendance_correlation_inverse_tendance devient sa propre decomposition de Cholesky L où L L^T est la matrice de depart.
  
//   gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,matrice_tendance_correlation_inverse_tendance,matrice_tendance_correlation_inverse); // matrice_tendance_correlation_inverse = matrice_tendance_correlation_inverse_tendance(deja choleskysee)^(-1) matrice_tendance_correlation_inverse

//   gsl_matrix* matrice_tendance_cholesky = gsl_matrix_alloc(tendance->size1, tendance->size2);
//   gsl_matrix_memcpy(matrice_tendance_cholesky,tendance); //matrice_tendance_cholesky = tendance
//   gsl_blas_dtrsm(CblasRight,CblasLower,CblasTrans,CblasNonUnit,1.0,matrice_tendance_correlation_inverse_tendance,matrice_tendance_cholesky); //matrice_tendance_cholesky = tendance matrice_tendnace_correlation_inverse_tendance(deja choleskysee)^(-T) 

//   gsl_matrix_free(matrice_tendance_correlation_inverse_tendance);

//   gsl_matrix* projectionBerger = gsl_matrix_alloc(matrice_correlation_derivee->size1, matrice_correlation_derivee->size2);
//   gsl_matrix_set_identity(projectionBerger); //projection = identite
//   gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1.0,matrice_tendance_cholesky,matrice_tendance_correlation_inverse,1.0,projectionBerger); //projectionBerger = projectionBerger(identite) - matrice_tendance_cholesky matrice_tendance_correlation_inverse

//   gsl_matrix_free(matrice_tendance_correlation_inverse);
//   gsl_matrix_free(matrice_tendance_cholesky);

//   gsl_blas_dsymm(CblasRight,CblasLower,1.0,matrice_correlation_derivee_copie,projectionBerger,0.0,matrice_correlation_derivee); // matrice_correlation_derivee = projectionBerger matrice_correlation_derivee


//         // FILE * h = fopen ("projectionBerger.txt", "w");
//         // gsl_matrix_fprintf (h, projectionBerger,"%f");
//         // fclose (h);

  
//   gsl_matrix_free(matrice_correlation_derivee_copie);

//   gsl_blas_dgemv(CblasNoTrans,1.0,projectionBerger,observations_sauvegarde,0.0,observations); // Dans la vraisemblance integree, il convient de remplacer les observations par projectionBerger %*% observations (c'est en fait y - H beta_chapeau). Attention, comme projectionBerger depend de theta, il est necessaire de garder en memoire les observations non modifiees : observations_sauvegarde
  
//   gsl_matrix_free(projectionBerger);
// }
  


//fin de Matern.cpp

//provient de MaternAnisotropeGeometrique.h



class NoyauMaternAnisotropeGeometrique : public NoyauMatern
{
 public :

  // NoyauMaternAnisotropeGeometrique(double const& regularite); // devenu INUTILE car repose sur NoyauMatern::NoyauMatern(double const& regularite)
  NoyauMaternAnisotropeGeometrique(double const& regularite, gsl_vector* longueurs_correlation);
  NoyauMaternAnisotropeGeometrique(NoyauMaternAnisotropeGeometrique const &noyau_a_copier);

  ~NoyauMaternAnisotropeGeometrique();

  double calCorrelation(gsl_vector* ecart_rd) const; // calcul de la correlation (suivant sa conception anisotrope geometrique) : necessite la donnee d'un ecart r-dimensionnel et de r longueurs de correlation. Le premier est donne en parametre, le second est membre de la classe mere.
  double calCorrelationDerivee(gsl_vector* ecart_rd, unsigned int const& indice_derive); // calcul de la derivee de la correlation par rapport a la indice_derive-ieme longueur de correlation.

};

// Sampling d'une longueur de correlation sachant toutes les autres par algorithme de Metropolis.
unsigned int choisitLongueurCorrelationMetropolis(const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, unsigned int const& indice_longueur_correlation_a_choisir, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);

// Sampling d'un nouveau jeu de longueurs de correlation par Gibbs.
  void generePointPosteriorGibbs(gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, double const& sd_instrum, unsigned int const& iterations, gsl_vector* vect_taux_sauts, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);

// Genere nb_points_a_generer jeux de longueurs de correlation par Gibbs.
void generePointsPosterior(RcppGSL::vector<double> & tous_les_points_posterior,unsigned int const& nb_points_a_generer, gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);



//fin de MaternAnisotropeGeometrique.h

//provient de MaternAnisotropeGeometrique.cpp



/*
NoyauMaternAnisotropeGeometrique::NoyauMaternAnisotropeGeometrique(double const& regularite) : NoyauMatern(regularite)
{
}
*/

NoyauMaternAnisotropeGeometrique::NoyauMaternAnisotropeGeometrique(double const& regularite, gsl_vector* longueurs_correlation) : NoyauMatern(regularite , longueurs_correlation)
{
}

NoyauMaternAnisotropeGeometrique::NoyauMaternAnisotropeGeometrique(NoyauMaternAnisotropeGeometrique const& noyau_a_copier) : NoyauMatern(noyau_a_copier)
{
}

NoyauMaternAnisotropeGeometrique::~NoyauMaternAnisotropeGeometrique()
{
  // RIEN DU TOUT : le destructeur de NoyauMatern (classe mere) s'occupe deja de liberer la memoire
}

double NoyauMaternAnisotropeGeometrique::calCorrelation(gsl_vector* ecart_rd) const
{
  if(getDimension()==1) // cas isotrope
  {
    gsl_vector_scale(ecart_rd,1.0/ getUneLongueurCorrelation(0)); // ecart_rd = ecart_rd / m_longueurs_correlation
  }
  
  else
  {
   gsl_vector_div(ecart_rd,m_longueurs_correlation); // ecart_rd = ecart_rd / m_longueurs_correlation
  }

  // Noyau anisotrope geometrique : calcule la correlation correspondant a la distance euclidienne modulee par les r longueurs de correlation correspondant aux r dimensions (rappel : r n'a pas ete formellement defini dans le code, il s'agit en fait de la taille de m_longueurs_correlation)
  double covar(1.0);
//  double distance(0);
  double distance(gsl_blas_dnrm2(ecart_rd));
//   for(int i=0; i<getDimension() ; i++) // Somme des carres des ecarts selon chaque dimension divises par la longueur de correlation correspondante.
//     {
//       distance += gsl_pow_2(gsl_vector_get(ecart_rd,i));
//     }
//   distance = sqrt(distance); // distance euclidienne - apres deformation de l'espace par les longueurs de correlation
  //cout << "LC_matcorrel = " << getUneLongueurCorrelation(0);
  //cout << "distance_matcorrel = " << distance;
  covar *= calCorrelation1d(distance);
  //cout << "covar_matcorrel = " << covar<<endl;
  return covar;
}

double NoyauMaternAnisotropeGeometrique::calCorrelationDerivee(gsl_vector* ecart_rd, unsigned int const& indice_derive)
{
  if(indice_derive>=getDimension())
    {
      cout << "Erreur dans NoyauMaternAnisotropeGeometrique::calCorrelationAnisGeomDerivee :" << endl;
      cout << "Indice derive : " << indice_derive << endl;
      cout << "Dimension     : " << getDimension() << endl;
      cout << "NB : Les indice_derive envisageables sont compris entre 0 et " << getDimension()-1 << endl;
      cout << "Je renvoie 0." << endl;
      return 0.0;
    }

  if(gsl_vector_get(ecart_rd,indice_derive)==0.0) // Si l'ecart selon la dimension qu'on derive est nul, la derivee est nulle.
    {
      return 0.0;
    }

//  cout<<"Taille de ecart_rd =" << ecart_rd->size;

// cout<<"ecart_rd =" << endl;
// for(int i=0; i<ecart_rd->size; i++)
// {
//   cout<<gsl_vector_get(ecart_rd,i)<< endl;
// }

//cout <<endl;

  if(getDimension()==1) // cas isotrope
  {
    gsl_vector_scale(ecart_rd,1.0/ getUneLongueurCorrelation(0)); // ecart_rd = ecart_rd / m_longueurs_correlation
    //cout<<"longcorr = " << getUneLongueurCorrelation(0);
  }
  else
  {
   gsl_vector_div(ecart_rd,m_longueurs_correlation); // ecart_rd = ecart_rd / m_longueurs_correlation
  }

// cout<<"ecart_rd apres division =" << endl;
// for(int i=0; i<ecart_rd->size; i++)
// {
//   cout<<gsl_vector_get(ecart_rd,i)<< endl;
// }

//  double distance(0);
  double distance(gsl_blas_dnrm2(ecart_rd));
//   for(int i=0; i<getDimension() ; i++) // Somme des carres des ecarts selon chaque dimension divises par la longueur de correlation correspondante.
//     {
//       distance += gsl_pow_2(gsl_vector_get(ecart_rd,i));
//     }
//   distance = sqrt(distance); // distance euclidienne - apres deformation de l'espace par les longueurs de correlation

  //cout << "distance = "<<distance;
  double covar(1.0);

  // La derivation par rapport a l'une des longueurs de correlation est compliquee dans le cas anisotrope geometrique. On n'a d'autre recours que d'ecrire la formule. Je ne vois pas de decomposition possible (contrairement a ce qui se passe dans le cas tensorise).
  double cste(2.0*sqrt(m_regularite));
  double num_standard_global(cste*distance); 
  if(num_standard_global>705.0) return 0.0;
  covar *= gsl_sf_gammainv(m_regularite) 
/  exp((m_regularite -1)* M_LN2)
* cste;

if(getDimension()==1) // cas isotrope
{
  covar *= distance;
  //cout<<"distance = " << distance;
  //cout<<"isotrope covar = "<<covar;
}
else 
{
  covar *= gsl_pow_2(gsl_vector_get(ecart_rd,indice_derive)) 
  / distance;
  //  cout<<"distance = " << distance;
  //cout<<"anisotrope covar = "<<covar;
}

covar *= 1.0 /gsl_vector_get(m_longueurs_correlation,indice_derive)
* exp( m_regularite * gsl_sf_log(num_standard_global) ) 
* gsl_sf_bessel_Knu(fabs(m_regularite-1), num_standard_global) ; // La fonction de Bessel modifiee de 2e espece K(nu,x) est telle que K(-nu,x)=K(nu,x).

//cout << "covar = " << covar;
  return covar;
}






unsigned int choisitLongueurCorrelationMetropolis(const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, unsigned int const& indice_longueur_correlation_a_choisir, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{
  double candidat;
  double densite_candidat;
  double rapport;
  unsigned int accepte;
  double init(1.0);

  unsigned int compte_sauts_acceptes(0);

    //La longueur de correlation a definir est initialisee a la valeur init.
  noyau->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,init);

  //Creation du noyau sur lequel on va travailler.
   NoyauMaternAnisotropeGeometrique* noyau_candidat = new NoyauMaternAnisotropeGeometrique(*noyau);
 

  // densite_init est la densite correspondant aux parametres initiaux
   double densite_init = densitePosterior(observations, indice_longueur_correlation_a_choisir, matrice_ecarts, noyau, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);


  // Sampling de la nouvelle longueur de correlation par Metropolis avec iterations pas.
  for(int i=0; i<iterations; i++)
    {
      candidat = noyau->getUneLongueurCorrelation(indice_longueur_correlation_a_choisir) + gsl_ran_gaussian(r,sd_instrum); // r est un gsl_rng*, objet mysterieux qui a un rapport avec le generateur aleatoire de nombres.

      //cout<<"Le candidat choisi est : "<<candidat<<endl;

      if(candidat>0.0)
  {
    noyau_candidat->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,candidat);

    //cout<<"Les longueurs_correlation sont " << noyau_candidat->getUneLongueurCorrelation(0) << " et " << noyau_candidat->getUneLongueurCorrelation(1) << endl << endl << endl;

    
    densite_candidat = densitePosterior(observations, indice_longueur_correlation_a_choisir, matrice_ecarts, noyau_candidat, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
    //cout<<"densite_candidat = "<<densite_candidat;
    //cout<<"densite_init = "<<densite_init;
    rapport = densite_candidat / densite_init;
    //cout<<"rapport = "<<rapport<<endl;
    accepte = gsl_ran_bernoulli(r,min(rapport,1.0));
    //cout<<"ACCEPTE ? "<<accepte<<endl;
    if(accepte==1)
      {
        //cout<<"rapport = "<<rapport<<endl; 
        //cout<<"ancienne valeur = " << noyau->getUneLongueurCorrelation(indice_longueur_correlation_a_choisir)<<endl;
        //cout << "candidat accepte = " << candidat << endl << endl;
    if(candidat<0.002) cout << "TRES PETIT !!!" << endl;
        //cout<<"On est entre dans accepte==1."<<endl;
        noyau->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,candidat);
        densite_init = densite_candidat;
        compte_sauts_acceptes++;
      }
  }
    }

  delete noyau_candidat;

  return compte_sauts_acceptes;
}



void generePointPosteriorGibbs(gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, double const& sd_instrum, unsigned int const& iterations, gsl_vector* vect_taux_sauts, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{
  double nouvelle_longueur_correlation;
  unsigned int nombre_sauts_acceptes;
  double iterations_double(iterations);

  double taux;

  // ATTENTION : IL FAUT FAIRE DU VRAI GIBBS, ie AVEC BALAYAGE ALEATOIRE

  double dbl_nombre_dimension = (double)point_posterior->size;
  double dbl_indice_longueur_correlation_a_choisir;

  int indice_longueur_correlation_a_choisir;
  int ancien_indice_longueur_correlation_a_choisir(-1);


  //LA BOUCLE EST FAUSSE
  // FAUX Pour chaque dimension, on choisit une nouvelle longueur de correlation.
  // FAUX for(int indice_longueur_correlation_a_choisir=0; indice_longueur_correlation_a_choisir<point_posterior->size; indice_longueur_correlation_a_choisir++)
  // FAUX  {

  // ON PEUT QUAND MEME FAIRE UNE BOUCLE : chaque nouveau sampling ne met a jour qu'une composante, donc il est raisonnable d'attendre quelques samplings avant d'enregistrer les changements dans le fichier de resultats
  for(int i=0; i<point_posterior->size; i++)
    {
dbl_indice_longueur_correlation_a_choisir =  gsl_ran_flat(r,0.0,dbl_nombre_dimension);
 indice_longueur_correlation_a_choisir = floor(dbl_indice_longueur_correlation_a_choisir);

 if(indice_longueur_correlation_a_choisir != ancien_indice_longueur_correlation_a_choisir) // Il est inutile de sampler plusieurs fois de suite le meme parametre
   {


      nombre_sauts_acceptes = choisitLongueurCorrelationMetropolis(r,observations, matrice_ecarts, noyau, indice_longueur_correlation_a_choisir, sd_instrum, iterations, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
taux = nombre_sauts_acceptes / iterations_double;

      
  gsl_vector_set(vect_taux_sauts, indice_longueur_correlation_a_choisir, taux);

      //cout<<"Marqueur2"<<endl;
      nouvelle_longueur_correlation = noyau->getUneLongueurCorrelation(indice_longueur_correlation_a_choisir);
      
      gsl_vector_set(point_posterior, indice_longueur_correlation_a_choisir, nouvelle_longueur_correlation);

      ancien_indice_longueur_correlation_a_choisir = indice_longueur_correlation_a_choisir;
   }
 
   }
}


void generePointsPosterior(RcppGSL::vector<double> & tous_les_points_posterior,unsigned int const& nb_points_a_generer, gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternAnisotropeGeometrique* noyau, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{

  gsl_vector* vect_taux_sauts = gsl_vector_alloc(point_posterior->size);

  for(int i=0; i<nb_points_a_generer; i++)
    {
      generePointPosteriorGibbs(point_posterior, r, observations, matrice_ecarts, noyau, sd_instrum, iterations, vect_taux_sauts, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);

      for(int j=0; j<point_posterior->size; j++) // remplit le vecteur tous_les_points_posterior
      {
        tous_les_points_posterior[i * (point_posterior->size) + j] = gsl_vector_get(point_posterior,j);
      }
   //       {  
   //   FILE * f = fopen ("point_posterior_tens.txt", "a");
   //   gsl_vector_fprintf (f, point_posterior,"%f");
   //   fclose (f);
   // }

   //       {  
   //   FILE * g = fopen ("taux_sauts.txt", "a");
   //   gsl_vector_fprintf (g, vect_taux_sauts,"%f");
   //   fclose (g);
   // }
    }    

  gsl_vector_free(vect_taux_sauts);
}
  
// fin de MaternAnisotropieGeometrique.cpp

//provient de MaternTensorise.h



class NoyauMaternTensorise : public NoyauMatern
{
 public :

  // NoyauMaternTensorise(double const& regularite); // devenu INUTILE car repose sur NoyauMatern::NoyauMatern(double const& regularite)
  NoyauMaternTensorise(double const& regularite, gsl_vector* longueurs_correlation);
  NoyauMaternTensorise(NoyauMaternTensorise const &noyau_a_copier);

  ~NoyauMaternTensorise();

  double calCorrelation(gsl_vector* ecart_rd) const; // calcul de la correlation (suivant sa conception tensorisee) : necessite la donnee d'un ecart r-dimensionnel et de r longueurs de correlation. Le premier est donne en parametre, le second est membre de la classe mere.
  double calCorrelationDerivee(gsl_vector* ecart_rd, unsigned int const& indice_derive); // calcul de la derivee de la correlation par rapport a la indice_derive-ieme longueur de correlation.

};

// Sampling d'une longueur de correlation sachant toutes les autres par algorithme de Metropolis.
unsigned int choisitLongueurCorrelationMetropolis(const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, unsigned int const& indice_longueur_correlation_a_choisir, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);

// Sampling d'un nouveau jeu de longueurs de correlation par Gibbs.
void generePointPosteriorGibbs(gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, double const& sd_instrum, unsigned int const& iterations, gsl_vector* vect_taux_sauts, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);

// Genere nb_points_a_generer jeux de longueurs de correlation par Gibbs.
void generePointsPosterior(RcppGSL::vector<double> & tous_les_points_posterior,unsigned int const& nb_points_a_generer, gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives); //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde);



//fin de MaternTensorise.h

//provient de MaternTensorise.cpp



/*
NoyauMaternTensorise::NoyauMaternTensorise(double const& regularite) : NoyauMatern(regularite)
{
}
*/

NoyauMaternTensorise::NoyauMaternTensorise(double const& regularite, gsl_vector* longueurs_correlation) : NoyauMatern(regularite, longueurs_correlation)
{
}


NoyauMaternTensorise::NoyauMaternTensorise(NoyauMaternTensorise const& noyau_a_copier) : NoyauMatern(noyau_a_copier)
{
}


NoyauMaternTensorise::~NoyauMaternTensorise()
{
  // RIEN DU TOUT : le destructeur de NoyauMatern (classe mere) s'occupe deja de liberer la memoire
}

double NoyauMaternTensorise::calCorrelation(gsl_vector* ecart_rd) const
{
  gsl_vector_div(ecart_rd, m_longueurs_correlation); // ecart_rd = ecart_rd / m_longueurs_correlation

  // Noyau tensorise : calcule le produit des correlations selon les r dimensions (rappel : r n'a pas ete formellement defini dans le code, il s'agit en fait de la taille de m_longueurs_correlation)
  double covar(1.0);
  for(int i=0; i<getDimension() ; i++)
    {
      covar *= calCorrelation1d(gsl_vector_get(ecart_rd,i));
      //cout<<"covar = "<<covar<<endl;
    }
  return covar;
}

double NoyauMaternTensorise::calCorrelationDerivee(gsl_vector* ecart_rd, unsigned int const& indice_derive)
{
  if(indice_derive>=getDimension())
    {
      cout << "Erreur dans NoyauMaternTensorise::calCorrelationDerivee :" << endl;
      cout << "Indice derive : " << indice_derive << endl;
      cout << "Dimension     : " << getDimension() << endl;
      cout << "NB : Les indice_derive envisageables sont compris entre 0 et " << getDimension()-1 << endl;
      cout << "Je renvoie 0." << endl;
      return 0.0;
    }

  if(gsl_vector_get(ecart_rd,indice_derive)==0.0) // Si l'ecart selon la dimension qu'on derive est nul, la derivee est nulle.
    {
      return 0.0;
    }

  gsl_vector_div(ecart_rd,m_longueurs_correlation); // ecart_rd = ecart_rd / m_longueurs_correlation


  // La derivee du produit est la derivee de la longueur de correlation correspondant a la dimension indice_derive multipliee par le produit des longueurs de correlation correspondant aux autres dimensions.
  double covar(1.0);

  // D'abord, calcul de la derivee de la longueur de correlation correspondant a la dimension indice_derive.
  double cste(2.0*sqrt(m_regularite));
  double num_standard_indice_derive(cste*fabs(gsl_vector_get(ecart_rd,indice_derive)));
  //cout << "num_standard_indice_derive = " << num_standard_indice_derive << endl;
  if(num_standard_indice_derive>705.0) return 0.0;
  covar *= gsl_sf_gammainv(m_regularite) 
/ gsl_sf_exp((m_regularite -1)* M_LN2)
/ gsl_vector_get(m_longueurs_correlation,indice_derive) 
* gsl_sf_exp( (m_regularite+1) * gsl_sf_log(num_standard_indice_derive) ) 
* gsl_sf_bessel_Knu(fabs(m_regularite-1), num_standard_indice_derive) ; // La fonction de Bessel modifiee de 2e espece K(nu,x) est telle que K(-nu,x)=K(nu,x).

  // cout<<"covar = "<<covar<<endl;

  // Puis, calcul du produit des longueurs de correlation correspondant aux autres dimensions.
  for (int i=1; i<getDimension(); i++) // Initialisation a 1 car indice_derive deja traite
    {
      covar *= calCorrelation1d(gsl_vector_get(ecart_rd,(indice_derive+i)%getDimension())); // on passe en revue les autres correlations, en commencant par la indice_derive+1 -ieme, et en terminant par la indice_derive-1 -ieme.
      // cout<<"covar = "<<covar<<endl;
    }

  return covar;
}





unsigned int choisitLongueurCorrelationMetropolis(const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, unsigned int const& indice_longueur_correlation_a_choisir, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{
  double candidat;
  double densite_candidat;
  double rapport;
  unsigned int accepte;
  double init(1.0);

  unsigned int compte_sauts_acceptes(0);


  //La longueur de correlation a definir est initialisee a la valeur init.
  noyau->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,init);

  //Creation du noyau sur lequel on va travailler.
    NoyauMaternTensorise* noyau_candidat = new NoyauMaternTensorise(*noyau);
  
  // densite_init est la densite correspondant aux parametres initiaux
    double densite_init = densitePosterior(observations, indice_longueur_correlation_a_choisir, matrice_ecarts, noyau, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);



  // Sampling de la nouvelle longueur de correlation par Metropolis avec iterations pas.
  for(int i=0; i<iterations; i++)
    {
      candidat = noyau->getUneLongueurCorrelation(indice_longueur_correlation_a_choisir) + gsl_ran_gaussian(r,sd_instrum); // r est un gsl_rng*, objet mysterieux qui a un rapport avec le generateur aleatoire de nombres.

      //cout<<"Le candidat choisi est : "<<candidat<<endl;

      //if(candidat <= 0.0) cout <<0 << endl;
      if(candidat>0.0)
  {
    noyau_candidat->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,candidat);

    //cout<<"Les longueurs_correlation sont " << noyau_candidat->getUneLongueurCorrelation(0) << " et " << noyau_candidat->getUneLongueurCorrelation(1) << endl << endl << endl;

    
    densite_candidat = densitePosterior(observations, indice_longueur_correlation_a_choisir, matrice_ecarts, noyau_candidat, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
    //cout<<"densite_candidat = "<<densite_candidat;
    //cout<<"densite_init = "<<densite_init;
    rapport = densite_candidat / densite_init;
    //cout<<"candidat = " << candidat << endl;
    //cout<<"rapport = "<<rapport<<endl;
    accepte = gsl_ran_bernoulli(r,min(rapport,1.0));
    //cout<<"ACCEPTE ? "<<accepte<<endl;
    if(accepte==1)
      {
        //cout<<"On est entre dans accepte==1."<<endl;
        noyau->modifieUneLongueurCorrelation(indice_longueur_correlation_a_choisir,candidat);
        densite_init = densite_candidat;

        compte_sauts_acceptes++;
      }

    //cout << rapport<<endl;
  }
    }
  //cout << nb_accepte << " deplacements acceptes." << endl;
  //if(nb_accepte==0) cout <<"Aucun deplacement accepte." << endl;
  //cout << "------------------------ Fin de la boucle for. -------------------" << endl;

  delete noyau_candidat;

  return compte_sauts_acceptes;
}



void generePointPosteriorGibbs(gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, double const& sd_instrum, unsigned int const& iterations, gsl_vector* vect_taux_sauts, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{
  double nouvelle_longueur_correlation;
  unsigned int nombre_sauts_acceptes;
  double iterations_double(iterations);

  double taux;

  // ATTENTION : IL FAUT FAIRE DU VRAI GIBBS, ie AVEC BALAYAGE ALEATOIRE

  double dbl_nombre_dimension = (double)point_posterior->size;
  double dbl_indice_longueur_correlation_a_choisir;

  int indice_longueur_correlation_a_choisir;
  int ancien_indice_longueur_correlation_a_choisir(-1);


  //LA BOUCLE EST FAUSSE
  // FAUX Pour chaque dimension, on choisit une nouvelle longueur de correlation.
  // FAUX for(int indice_longueur_correlation_a_choisir=0; indice_longueur_correlation_a_choisir<point_posterior->size; indice_longueur_correlation_a_choisir++)
  // FAUX  {

  // ON PEUT QUAND MEME FAIRE UNE BOUCLE : chaque nouveau sampling ne met a jour qu'une composante, donc il est raisonnable d'attendre quelques samplings avant d'enregistrer les changements dans le fichier de resultats
  for(int i=0; i<point_posterior->size; i++)
    {
dbl_indice_longueur_correlation_a_choisir =  gsl_ran_flat(r,0.0,dbl_nombre_dimension);
 indice_longueur_correlation_a_choisir = floor(dbl_indice_longueur_correlation_a_choisir);

 if(indice_longueur_correlation_a_choisir != ancien_indice_longueur_correlation_a_choisir) // Il est inutile de sampler plusieurs fois de suite le meme parametre
   {


      nombre_sauts_acceptes = choisitLongueurCorrelationMetropolis(r,observations, matrice_ecarts, noyau, indice_longueur_correlation_a_choisir, sd_instrum, iterations, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
taux = nombre_sauts_acceptes / iterations_double;

      
  gsl_vector_set(vect_taux_sauts, indice_longueur_correlation_a_choisir, taux);

      //cout<<"Marqueur2"<<endl;
      nouvelle_longueur_correlation = noyau->getUneLongueurCorrelation(indice_longueur_correlation_a_choisir);
      
      gsl_vector_set(point_posterior, indice_longueur_correlation_a_choisir, nouvelle_longueur_correlation);

      ancien_indice_longueur_correlation_a_choisir = indice_longueur_correlation_a_choisir;
   }
 
   }
}



void generePointsPosterior(RcppGSL::vector<double> & tous_les_points_posterior, unsigned int const& nb_points_a_generer, gsl_vector* point_posterior, const gsl_rng* r, gsl_vector* observations, MatriceEcarts* matrice_ecarts, NoyauMaternTensorise* noyau, double const& sd_instrum, unsigned int const& iterations, const gsl_matrix* injection_orthogonal_tendance, int const& nombre_observations_effectives) //bool const& berger, const gsl_matrix* tendance, const gsl_vector* observations_sauvegarde)
{

  gsl_vector* vect_taux_sauts = gsl_vector_alloc(point_posterior->size);

  for(int i=0; i<nb_points_a_generer; i++)
    {
      generePointPosteriorGibbs(point_posterior, r, observations, matrice_ecarts, noyau, sd_instrum, iterations, vect_taux_sauts, injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);

      for(int j=0; j<point_posterior->size; j++) // remplit le vecteur tous_les_points_posterior
      {
        tous_les_points_posterior[i * (point_posterior->size) + j] = gsl_vector_get(point_posterior,j);
      }
   //       {  
   //   FILE * f = fopen ("point_posterior_tens.txt", "a");
   //   gsl_vector_fprintf (f, point_posterior,"%f");
   //   fclose (f);
   // }

   //       {  
   //   FILE * g = fopen ("taux_sauts.txt", "a");
   //   gsl_vector_fprintf (g, vect_taux_sauts,"%f");
   //   fclose (g);
   // }
    }

  gsl_vector_free(vect_taux_sauts);
}
  
// provient de main.cpp



#define fichier_observations "observations.txt"
#define fichier_planXP "planXP.txt"
#define fichier_tendance "tendance.txt"
#define fichier_metaparametres "metaparametres.txt"

#define fichier_posterior_anis_geom "point_posterior_anis_geom.txt"
#define fichier_posterior_tens "point_posterior_tens.txt"

#define fichier_type_noyau "type_noyau.txt"
#define fichier_type_prior "type_prior.txt"

gsl_rng * r;  /* global generator */

// Nombre de points du plan d'experience (planXP) devine grace au nombre de lignes du fichier qui les recense.
int compteLignes(string nom_fichier)
{
ifstream f(nom_fichier.c_str());
string line;
int compteur_lignes;
for (compteur_lignes = 0; getline(f, line); ++compteur_lignes) ;
f.close();
return compteur_lignes;
}


int compteColonnes(string nom_fichier)
    {
    ifstream myfile(nom_fichier.c_str());
    string line, temp;
    stringstream ss;
    int ncols=0;
    
    //discard first line from the file (it's a header)
    // NON, ce n'en est pas un
    //getline(myfile, line);
    


    getline(myfile, line);
    ss.clear();
    ss << line;
    
while (ss >> temp)
    {
        ncols++;
    }

return ncols;
    }



//[[Rcpp::export]]
RcppGSL::vector<double> GibbsPosteriorC(int dimension, double regularite, int nb_points_a_generer, int iterations, double sd_instrum, int nb_fonctionsTendance, unsigned long int random_seed, std::vector<std::string> words, RcppGSL::matrix<double> planXP, RcppGSL::vector<double> observations, RcppGSL::matrix<double> tendanceR)
{
 
  // Je ne comprends pas grand chose a ce qui suit : j'ai juste suivi le modele donne par la documentation de gsl_randist.
  //********* Preparation du generateur aleatoire de nombres**************
  // const gsl_rng_type * T;
  // gsl_rng_env_setup();
  // T = gsl_rng_default;
  // r = gsl_rng_alloc (T);
  //****** Fin de la preparation du generateur aleatoire de nombres******

// *********Une maniere plus intelligente d'initialiser le generateur aleatoire
  r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, random_seed);
  //****** Fin de la preparation intelligente du generateur aleatoire de nombres******

  // Les valeurs observees sont lues et stockees dans le gsl_vector* observations.
//   FILE * fObservations = fopen(fichier_observations,"r");
// gsl_vector* observations = gsl_vector_alloc(compteLignes(fichier_observations));
//   if(fObservations != NULL)
//       {
//   gsl_vector_fscanf(fObservations,observations);
//   fclose(fObservations);
//       }


  
  // L'utilisateur renseigne les metaparametres dans fichier_metaparametres. Ces metaparametres sont recuperes d'abord dans un gsl_vector*, puis separement dans des int et double.
// FILE * fMetaparametres = fopen(fichier_metaparametres,"r");
// gsl_vector* metaparametres = gsl_vector_alloc(6);
// gsl_vector_fscanf(fMetaparametres,metaparametres);
// double regularite = gsl_vector_get(metaparametres,1);
// int dimension = floor(gsl_vector_get(metaparametres,0));
// int nb_points_a_generer = floor(gsl_vector_get(metaparametres,2));
//  int iterations = floor(gsl_vector_get(metaparametres,3));
//  double sd_instrum = gsl_vector_get(metaparametres,4);
// int nb_fonctionsTendance = gsl_vector_get(metaparametres,5);
// gsl_vector_free(metaparametres);
// fclose(fMetaparametres);

RcppGSL::vector<double> tous_les_points_posterior(dimension * nb_points_a_generer);

// Les emplacements des points du planXP sont lus dans un fichier presente comme une matrice (lecture d'abord de gauche a droite puis de haut en bas) et stockes dans la gsl_matrix* planXP.
// FILE * fPlanXP = fopen(fichier_planXP,"r");
// gsl_matrix* planXP = gsl_matrix_alloc(compteLignes(fichier_planXP),compteColonnes(fichier_planXP));      //,dimension);
// gsl_matrix_fscanf(fPlanXP,planXP);

// Allocation du gsl_vector* tampon point_posterior qui stockera le point (ie le vecteur contenant les longueurs de correlation) en train d'être sample a partir de la loi a posteriori jointe.
gsl_vector* point_posterior = gsl_vector_alloc(dimension);

// Allocation et remplissage de la MatriceEcarts (classe contenant... la matrice des ecarts entre points du planXP)
MatriceEcarts* matrice_ecarts = new MatriceEcarts(planXP);
cout << "Dimension de matrice_ecarts = " << matrice_ecarts->getDimension()<<endl;

// Allocation d'un vecteur contenant les valeurs initiales des longueurs de correlation.
gsl_vector* longueurs_correlation = gsl_vector_alloc(dimension);
gsl_vector_set_all(longueurs_correlation,1.0);


// Attention !!! On peut avoir un probleme si, dans la partie en R, tendance n'a pas ete defini comme une matrice (cas du krigeage simple)

// Les fonctions de base de la tendance sont lues dans un fichier presente comme une matrice (lecture d'abord de gauche a droite puis de haut en bas) et stockes dans la gsl_matrix* tendance.
// gsl_matrix * tendance = NULL;
// if(nb_fonctionsTendance > 0)
// {
// FILE * fTendance = fopen(fichier_tendance,"r");
// tendance = gsl_matrix_alloc(compteLignes(fichier_tendance),nb_fonctionsTendance);
// gsl_matrix_fscanf(fTendance,tendance);
// }

gsl_matrix* tendance = NULL;

if(tendanceR->size1 > 1) // Si tendance a une seule ligne, cela signifie ou bien que planXP a un seul point (impossible) ou bien que tendance est la matrice 1x1 nulle, ce qui code l'absence de fonction de tendance (krigeage simple)
{
  tendance = tendanceR;
}

cout << tendance << endl;

// Noyau de Matern anisotrope geometrique ou tensorise ? Reponse dans fichier_type_noyau.
//    vector<string> words; // peut devenir utile s'il est un jour necessaire de lire plusieurs mots dans le fichier
//ifstream fType(fichier_type_noyau);
//  string word;
//  while(fType >> word) { // boucle necessaire s'il faut lire plusieurs mots
//    words.push_back(word); // s'il faut lire plusieurs mots
//    } // s'il faut lire plusieurs mots
//  fType.close();
    string type = words[0]; // s'il faut lire plusieurs mots
  string type_prior("par_defaut_entraine_fin_du_programme");
  if(words.size()>1) type_prior = words[1];


//cout << "words = " << endl;
//cout << "type_prior = " << type_prior << endl;

 string nom_fichier;
 if(type == "geometrique") nom_fichier = fichier_posterior_anis_geom;
 else if (type == "tensorise") nom_fichier = fichier_posterior_tens;
 else
   {
     cout << "Je ne comprends pas l'instruction : geometrique ou tensorise ? Fin du programme." << endl;
     return 1;
   }


bool berger;
 int nombre_observations_effectives(0);
if(type_prior == "berger") 
{
  berger = 1;
  cout << "Le prior de Berger est utilise." << endl;
  nombre_observations_effectives = observations->size;
}
else if(type_prior == "REML") 
{
  berger = 0;
  cout << "Le prior utilisant la REstricted Likelihood est utilise." << endl;
  nombre_observations_effectives = observations->size - nb_fonctionsTendance;
}
else
   {
     cout << "Je ne comprends pas l'instruction : berger ou REML ? Fin du programme." << endl;
     return 1;
   }

// if(berger & nb_fonctionsTendance < 1)
// {
//  cout << "Avec le prior de Berger, il faut au moins une fonction de base pour la tendance." << endl;
//  return 1;
// }

// Definition du pointeur vers la matrice qui servira a la projection (prior berger ou REML)
 gsl_matrix* injection_orthogonal_tendance = NULL;
 // if(berger)
 //   {
 //     injection_orthogonal_tendance = gsl_matrix_alloc(planXP->size1,planXP->size1);
 //   }
 // else
 //   {
     injection_orthogonal_tendance =  gsl_matrix_alloc(planXP->size1,planXP->size1 - nb_fonctionsTendance);
 // }
 
     calInjectionREML(tendance,injection_orthogonal_tendance); //,berger);


 // Les observations doivent etre amputees de la partie qui sert a determiner la moyenne
 gsl_vector* observations_orthogonal_tendance = gsl_vector_alloc(injection_orthogonal_tendance->size2);
 gsl_blas_dgemv(CblasTrans,1.0,injection_orthogonal_tendance,observations,0.0,observations_orthogonal_tendance);
 gsl_vector* observations_sauvegarde = observations;
 observations = observations_orthogonal_tendance;
 gsl_vector_free(observations_sauvegarde);


 // Attention !!! Commenter ce qui suit empeche de specifier ou commence le Gibbs : on a une perte de fonctionnalite


 // FILE * fLC;

 // fLC = fopen(nom_fichier.c_str() ,"r");
 //    if(fLC != NULL)
 //      {
 //  int nombre_lignes = compteLignes(nom_fichier);
 //  int nombre_points = nombre_lignes / dimension;
 //  if(nombre_points == 0)
 //    {
 //      fclose(fLC);
 //      remove(nom_fichier.c_str());
 //    }
  
 //  else
 //    {
 //      gsl_vector* toutes_les_coordonnees = gsl_vector_alloc(nombre_points * dimension );
 //      gsl_vector_fscanf(fLC,toutes_les_coordonnees);
 //      fclose(fLC);
 //      fLC =  fopen(nom_fichier.c_str(),"w+");
 //      gsl_vector_fprintf(fLC,toutes_les_coordonnees,"%f");
 //      fclose(fLC);
 //      cout << "Longueurs de correlation initiales :" << endl;
 //      for(int i=(nombre_points-1)*dimension; i< nombre_points * dimension; i++)
 //        {
 //    gsl_vector_set(longueurs_correlation,i%dimension,gsl_vector_get(toutes_les_coordonnees,i));
 //    cout << "Numero " << i%dimension +1 << " : " << gsl_vector_get(longueurs_correlation,i%dimension) << endl;
 //        }
 //      gsl_vector_free(toutes_les_coordonnees);
 //    }
 //      }

    //cout<<"dans main LC = " << gsl_vector_get(longueurs_correlation,0)<<endl;

    
// Creation du noyau de Matern - anisotrope geometrique ou tensorise - qui permettra de construire les matrices de correlation.
    if(type == "geometrique")
      {
  cout << "Noyau de Matern anisotrope geometrique de dimension " << dimension << "." << endl;
  NoyauMaternAnisotropeGeometrique* noyauG = new NoyauMaternAnisotropeGeometrique(regularite,longueurs_correlation);
  //cout << "dans main LC noyauG = " << noyauG->getUneLongueurCorrelation(0)<<endl;
  generePointsPosterior(tous_les_points_posterior,nb_points_a_generer,point_posterior,r,observations,matrice_ecarts,noyauG,sd_instrum,iterations,injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
  delete noyauG;
      }
    else if(type == "tensorise")
      {
  cout << "Noyau de Matern tensorise de dimension " << dimension << "." << endl;  
  NoyauMaternTensorise* noyauT = new NoyauMaternTensorise(regularite,longueurs_correlation);
  //cout << noyauT->calCorrelation1d(287.95) << endl;
  generePointsPosterior(tous_les_points_posterior,nb_points_a_generer,point_posterior,r,observations,matrice_ecarts,noyauT,sd_instrum,iterations,injection_orthogonal_tendance,nombre_observations_effectives); //berger,tendance,observations_sauvegarde);
  delete noyauT;
      }
    else
      {
  cout << "Je n'ai pas compris quel type de noyau doit etre utilise." << endl;
      }
    
 // Liberation de la memoire (inutile en fait vu qu'on arrive a la fin du programme).
// gsl_vector_free(observations); //  J'AI PEUR QUE CELA SUPPRIME y_connus dans la memoire de R
// gsl_vector_free(observations_sauvegarde);
//gsl_matrix_free(planXP); //  J'AI PEUR QUE CELA SUPPRIME x_connus dans la memoire de R
//gsl_matrix_free(tendance); // J'AI PEUR QUE CELA SUPPRIME tendance dans la memoire de R
gsl_vector_free(longueurs_correlation);
 gsl_vector_free(point_posterior);
delete matrice_ecarts;


  return tous_les_points_posterior;

}

//[[Rcpp::export]]
int quarantedeux()
{
  return 42;
}
