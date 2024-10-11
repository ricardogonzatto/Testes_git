#33include <stdio.h> OPA
#33include <stdlib.h>
#include <string.h> quem ta lendo é gay
#include <math.h>
#include <gmp.h>
#include <mps/mps.h>
#include <nlopt.h>
#include <complex.h>
#include <time.h>
#include <chealpix.h>

#define LMAX 2000
#define NSIDE 64
#define NUMERO_SIMULACOES 1202
#define MVS_NUMERO_LINHAS (((LMAX-1)+1) * ((LMAX-1)+2)) - 2


typedef struct _FrechetData
{
  int l;
  double * restrict x;
  double * restrict y;
  double * restrict z;
} FrechetData;

#define ACOS(ab) acos (((ab) > 1.0) ? 1.0 : (((ab) < -1.0) ? -1.0 : (ab)))

void 
guess(int ell, const double *x, const double *y, const double *z, double *s)
{
//    char psis_write[300] = "psis.dat";  //psis_file
    int npix = (12 * NSIDE * NSIDE) / 2; // Total number of pixels divide per 2 
//    int npix = (12 * NSIDE * NSIDE);

    /*All ipix that we work*/
    double pixel_coords[npix][3];
    double vec[3];

    for (int ipix = 0; ipix < npix; ipix++) {
        double vec[3];
        pix2vec_ring(NSIDE, ipix, vec);
        pixel_coords[ipix][0] = vec[0];
        pixel_coords[ipix][1] = vec[1];
        pixel_coords[ipix][2] = vec[2];
    }

    /*This for calculate the psi function of MVs for each ipix and return the "chutes"*/
    double psi_min = 1.0e300;
    double chute_x, chute_y, chute_z;
    double psi[npix];
    double guess[3];

    for(int ipix = 0; ipix < npix; ipix ++){
        double sum_arccos_squared = 0.0;
        double dot_product = 0.0;
        
        for(int pos_mv = 0; pos_mv < (ell); pos_mv ++){

            dot_product =  (pixel_coords[ipix][0] * x[pos_mv]) + (pixel_coords[ipix][1] * y[pos_mv]) + (pixel_coords[ipix][2] * z[pos_mv]);
        //    sum_arccos_squared += ACOS(dot_product) * ACOS(dot_product);
            sum_arccos_squared += (ACOS(dot_product) * ACOS(dot_product)) + ((M_PI - ACOS(dot_product)) * (M_PI - ACOS(dot_product)));
        }
        
        psi[ipix] = sum_arccos_squared;
 //       printf("ipix = %d and PSI = %f\n", ipix, sum_arccos_squared);
 //       printf("psi_min = %f\n", psi_min);
        if(sum_arccos_squared < psi_min){

            psi_min = sum_arccos_squared;
            
            guess[0] = pixel_coords[ipix][0];
            guess[1] = pixel_coords[ipix][1];
            guess[2] = pixel_coords[ipix][2];
        }
    
    }
//    printf("For ell = %d: Valores mínimos encontrados (%lf, %lf, %lf)\n", ell, guess[0], guess[1], guess[2]);
    double theta_frechet[1], phi_frechet[1];
    vec2ang(guess, theta_frechet, phi_frechet);

    s[0] = theta_frechet[0];
    s[1] = phi_frechet[0];
//    printf("For ell = %d: Valores mínimos encontrados (%lf, %lf)\n", ell, s[0], s[1]);
    //printf("For ell = %d: Valores mínimos encontrados (%lf, %lf, %lf)\n", ell, chute_x, chute_y, chute_z);

    /*Saving the PSIS in ipix order*/
/*    FILE *PSIs = fopen(psis_write, "a");
    for(int ipix = 0; ipix < npix; ipix ++){

        fprintf(PSIs, "%.15f\n", psi[ipix]);
    }
    fclose(PSIs); */

}


double 
frechet_pol_min (unsigned n, const double *x, double *grad, void *my_func_data)
{
  FrechetData *fd = (FrechetData *) my_func_data;
  const int twol = 2 * fd->l;
//  const int twol = fd->l;
  double fechet_mean = 0.0;  
  double x_v, y_v, z_v;
  int i;
  
  {
    const double theta = x[0];
    const double phi   = x[1];
    
    x_v = sin (theta) * cos (phi);
    y_v = sin (theta) * sin (phi);
    z_v = cos (theta);

    for (i = 0; i < twol; i++)
    {
      const double ab      = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
      const double acos_ab = ACOS (ab);
      fechet_mean += acos_ab * acos_ab;
    }
    //printf (".", fechet_mean);
    if (grad != NULL)
    {
      const double dxdtheta = cos (theta) * cos (phi);
      const double dydtheta = cos (theta) * sin (phi);
      const double dzdtheta = - sin (theta);
    
      const double dxdphi = - sin (theta) * sin (phi);
      const double dydphi = + sin (theta) * cos (phi);
    
      grad[0] = 0.0;
      grad[1] = 0.0;
    
      for (i = 0; i < twol; i++)
      {
        const double ab      = (x_v * fd->x[i]) + (y_v * fd->y[i]) + (z_v * fd->z[i]);
        const double acos_ab = ACOS (ab);
        const double sqrt_1mab = sqrt (1.0 - ab * ab);
      
        grad[0] -= 2.0 * acos_ab * (dxdtheta * fd->x[i] + dydtheta * fd->y[i] + dzdtheta * fd->z[i]) / sqrt_1mab;
        grad[1] -= 2.0 * acos_ab * (dxdphi   * fd->x[i] + dydphi   * fd->y[i]                      ) / sqrt_1mab;
      }
    }
  }

  return fechet_mean;
}


void
frechet_pol2 (int l, double * restrict theta, double * restrict phi, double *frechet_vec_theta, double *frechet_vec_phi)
{
//  double frechet_vec_theta;
//  double frechet_vec_phi;
  const int twol = 2 * l;
  double x[2 * l], y[2 * l], z[2 * l];
  double x_ant[2 * l], y_ant[2 * l], z_ant[2 * l];
  double f;
  int i;
  FrechetData fd = {l, x, y, z};
  
  *frechet_vec_theta = 0.0;
  *frechet_vec_phi = 0.0;


  double ant_theta[l], ant_phi[l];

  int muda = 0;
  for (i = 0; i < twol; i++)
  {
    if(theta[i] < M_PI/2)
    {
      ant_theta[muda] = theta[i];
      ant_phi[muda] = phi[i];
      muda++;
      
    }
  }
//  printf("final = %d, sendo ell = %d \n", muda, l);
  if(muda != l){
    printf("Antipodes Number not correspond\n");
  }

  for (i = 0; i < l; i++)
  {
    x_ant[i] = sin (ant_theta[i]) * cos (ant_phi[i]);
    y_ant[i] = sin (ant_theta[i]) * sin (ant_phi[i]);
    z_ant[i] = cos (ant_theta[i]);
  }


  for (i = 0; i < twol; i++)
  {
    x[i] = sin (theta[i]) * cos (phi[i]);
    y[i] = sin (theta[i]) * sin (phi[i]);
    z[i] = cos (theta[i]);
  }

  {
//    double lb[2] = {0,  0};
//    double ub[2] = {M_PI,  + 2.0 * M_PI};

//    double lb[2] = {-1.0 * M_PI,  -2.0 * M_PI};
//    double ub[2] = {+1.0 * M_PI,  +2.0 * M_PI};
    double lb[2] = {0,  0};
    double ub[2] = {M_PI,  + 2.0 * M_PI};
    double s[2]  = {0.0, 0.0};
    double min_f = 1.0e300;
//    int ntheta = 8;
//    int nphi   = 8;
//    int ntheta = 15;
//    int nphi   = 15;
    

    nlopt_opt opt;
    
    opt = nlopt_create (NLOPT_LN_NELDERMEAD, 2);

    nlopt_set_lower_bounds (opt, lb);
    nlopt_set_upper_bounds (opt, ub);
    nlopt_set_min_objective(opt, &frechet_pol_min, &fd);
    nlopt_set_xtol_rel(opt, 1.0e-7);
//    nlopt_set_xtol_rel(opt, 1.0e-1);

    guess(l, x_ant, y_ant, z_ant, s);

//    s[0] = 1.0 * M_PI / ((ntheta) + 1.0) * (i + 1.0);
//    s[1] = 2.0 * M_PI / ((nphi)   + 1.0) * (j + 1.0);

    if (nlopt_optimize(opt, s, &f) < 0) 
    {
      printf("nlopt failed!\n");
    }
    else 
    {
//          printf("M2 % 22.15g % 22.15g % 22.15g\n", s[0], s[1], f);
      if (f < min_f)
      {
        min_f = f;

          *frechet_vec_theta = s[0];
          *frechet_vec_phi   = s[1];
      }
    }


    if (*frechet_vec_theta > (M_PI / 2))
  {
    *frechet_vec_theta = M_PI - *frechet_vec_theta;
    *frechet_vec_phi   = M_PI + *frechet_vec_phi;

    if (*frechet_vec_phi > (2 * M_PI))
      *frechet_vec_phi = *frechet_vec_phi - (2 * M_PI);
  }

//   printf ("M2 % 22.15g % 22.15g % 22.15g\n", *frechet_vec_theta, *frechet_vec_phi, min_f);
    nlopt_destroy (opt);

//    *frechet_vec_eta = 1-cos(frechet_vec_theta); 
//    *frechet_vec_varphi = frechet_vec_phi/(2*M_PI);
  } 
}


int main (int argc, char *argv[ ])

{
    char mvs_read[300], fvs_write[300];
    double frechet_vec_theta, frechet_vec_phi;
    static double theta[MVS_NUMERO_LINHAS], phi[MVS_NUMERO_LINHAS];

    int mc = atoi(argv[1]);

    sprintf(mvs_read,"/share/storage5/ricardogonzatto/NEW_MVS/fullsky/mvs_1000mc_control/mvstheory_fullsky_ellmax2000_mc%d.dat", mc);
    sprintf(fvs_write,"/share/storage5/ricardogonzatto/NEW_FVS/1000mc_theoric/fullsky/fvscontrol_fullsky_ellmax2000_mc%d.dat", mc);
    
    FILE *MVs = fopen(mvs_read, "r");
    for (int j = 0; j < MVS_NUMERO_LINHAS; j++)
    {
        fscanf(MVs, "%lf %lf\n", &theta[j], &phi[j]);
      //printf("linha = %d\n", j);
    }
    fclose (MVs);

    FILE *FVs = fopen(fvs_write, "a");
    for (int l = 2; l <= (LMAX); l++ )
    {
        int lim_min = ((pow((l-1),2)) + (l-1) - 2) / 2 ;
        int lim_max = ((pow(l,2)) + l - 2) / 2;

        double theta_each_l[2*l], phi_each_l[2*l];  
        int n = 0;

        for (int pos = lim_min; pos < lim_max; pos++ )
        {
            theta_each_l[n] = theta[pos];
            phi_each_l[n] = phi[pos];
            
            
            //Aqui eu adiciono as antipodas

            if(theta[pos] >= 0 && theta[pos] <= M_PI/2)
            {
               theta_each_l[l + n] = M_PI - theta[pos];
            } 

            if(phi[pos] >= 0 && phi[pos] < 2*M_PI)
            {
              phi_each_l[l + n] = fmod((phi[pos] + M_PI), (2*M_PI));
            }

            if(phi[pos] >= -M_PI && phi[pos] < M_PI)
            {
              phi_each_l[l + n] = fmod(((phi[pos] + M_PI) + M_PI), (2*M_PI)) - M_PI;
            } 

            n++;
            //
        }
        
        frechet_pol2 (l, theta_each_l, phi_each_l,  &frechet_vec_theta, &frechet_vec_phi);
        //printf("Theta = %f  Phi = %f\n", frechet_vec_theta, frechet_vec_phi);
        fprintf (FVs, "%.15f %.15f\n", frechet_vec_theta, frechet_vec_phi);

    }
    fclose (FVs);

    


}
