#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pas 0.00001 //On change les occurence de "pas" par "0.00001".

double x2_x_1(double x) {
    return x*x - x -1;
}

double dichotomie(double f(double), double x, double y, int N){
    if (f(x)*f(y) > 0){ //On vérifié que notre intervalle de départ est correct
	printf("Wrong interval\n");
	return 0.;
    }
    double z;
    int i;
    for(i = 0; i<N; i++){
        z = (x + y) /2.; //On se place au milieu entre x et y
        if (f(z)*f(x)<0) y = z; // si f(x)f(z) est négatif, z est du signe de y
        else x = z; //sinon il est du signe de x
    }
    return z; //On pourrait aussi renvoyer x ou y ça change pas grand chose.
}

double deriv(double f(double), double x){
    return( (f(x+1e-10) - f(x)) /1e-10 ); // Ici epsilon vaut 10^-10. On pourrait aussi le mettre en paramètre.
}

double newton(double f(double), double x, int N){
    int i;
    for(i = 0; i<N; i++){
        x = -f(x)/deriv(f,x)+x; // solution de f'(xn) (x - xn) +f(xn) = 0
    }
    return x;

}

double rectangle(double f(double), double a, double b, int N){
    double res =0; //On initie l'intégrale à 0
    int i;
    for(i=0; i<N; i++){
        res += (b-a)/N * f(a + i*(b-a)/N); //On ajoute l'air du rectangle de largeur (b-a)/N et d'hauteur f(a +i*(b-a)/N)
    }
    return res;
}


double trapeze(double f(double), double a, double b, int N){
    double res =0;
    double larg = (b-a)/N;
    for(int i=0; i<N; i++){
        res += f(a+ i*larg) + f(a+ (i+1)*larg);
    }
    res = res*larg/2;
    return res;
}


double max(double f(double), double a, double b, double eps){
    double maxi = f(a);
    for (double x = a; x<b; x+= eps){ //Notez ici la formulation de la boucle
        if (f(x)>maxi) maxi = f(x);
    }
    return maxi;

}

double min(double f(double), double a, double b, double eps){
    double mini = f(a);
    for (double x = a; x<b; x += eps){
        if (f(x)<mini) mini = f(x);
    }
    return mini;
}

double monte_carlaux(double f(double), double a, double b, int N){ //Ne fonctionne que pour des fonctions positives
    int i;
    double maxf = max(f,a,b,pas); //On récupère le max de f
    double sous_la_courbe = 0; // Variable auxiliaire comptant le nombre de points sous la courbe
    for(i = 0; i<N; i++)
    {
        double x = (double) rand();
        double y = (double) rand(); //x et y sont des nombres aléatoires
        double posx = a + (b-a)* x/RAND_MAX; //On génère depuis x une position aléatoire entre a et b
        double posy = maxf * (double) y/RAND_MAX; //On génère depuis y une position aléatoire entre 0 et maxf (notez que c'est la même ligne que ci-dessus avec a =0)
        if (f(posx) > posy) sous_la_courbe +=1; //Si le point est sous la courbe on incrémente sous_la_courbe
    }
    return (double) sous_la_courbe/N*maxf*(b-a); //La proportion de point sous la courbe est sous_la_courbe/N. L'air du rectangle est maxf*(b-a)
}

double monte_carlo(double f(double), double a, double b, int N){
   double minf = min(f,a,b,pas); //On calcule le minimum de f
   double f_aux(double a){ //On définit une fonction auxilliaire toujours positive
	return f(a) - minf ;
	}
   double int_f_aux = monte_carlaux(f_aux,a,b,N); // On calcule l'intégrale de f_aux
   return (b-a)*minf + int_f_aux; //On renvoie l'intégrale de f par linéarité de l'intégrale
}

int main()
{
    double sol = monte_carlo(x2_x_1,0.,3.,10000000);
    printf("%.3f\n",sol); //"%.3f s'arrete aux 3 premiers chiffres après la virgule.
    return 0;
}
