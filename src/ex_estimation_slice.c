#include "util.h" 

void mcmc(double **data1,double **data2, int *n1, int *n2, int *nb_col1,
	  int *B, int *dye_swap,
	  double *gamma_1, double *gamma1_p, 
	  double *gamma_2, double *gamma2_p, 
	  double *mu, double *mu_p,
	  double *beta2, double *beta2_p,
	  double *alpha2, double *alpha2_p,
	  double *delta22,
	  double *delta22_p,
	  double *eta, double *eta_p,
	  double *lambda_eps1, double *lambda_eps1_p,
	  double *lambda_eps2, double *lambda_eps2_p,
	  double *a_eps, double *b_eps,
	  double *a_eps_p, double *b_eps_p,
	  double *w, double *df_choice, int *nb_df,
	  double *df, double *df_p, double *w_p,
	  double *rho, double *rho_p,
	  double *lambda_gamma1, double *lambda_gamma2,
	  double *lambda_gamma1_p, double *lambda_gamma2_p,
 	  int *min_iter, int *batch, int *all_out, int *verbose);

void ex_R_link_mcmc(double *data_vec1, double *data_vec2, 
		    int *n1, int *n2, int *nb_col1, int *B, int *dye_swap,
		    double *gamma_1, double *gamma1_p, 
		    double *gamma_2, double *gamma2_p, 
		    double *mu, double *mu_p,
		    double *beta2, double *beta2_p,
		    double *alpha2, double *alpha2_p,
		    double *delta22,
		    double *delta22_p,
		    double *eta, double *eta_p,
		    double *lambda_eps1, double *lambda_eps1_p,
		    double *lambda_eps2, double *lambda_eps2_p,
		    double *a_eps, double *b_eps,
		    double *a_eps_p, double *b_eps_p,
		    double *w, double *df_choice, int *nb_df,
		    double *df, double *df_p, double *w_p,
		    double *rho, double *rho_p,
		    double *lambda_gamma1, double *lambda_gamma2,
		    double *lambda_gamma1_p, double *lambda_gamma2_p,
		    int *min_iter, int *batch, int* all_out, int *verbose)
{
  
  double **data1;
  double **data2;
  
  GetRNGstate();
  data1=dmatrix(*n1, *n2);  
  vec_mat(data_vec1,n1,n2,data1);
  data2=dmatrix(*n1, *n2);  
  vec_mat(data_vec2,n1,n2,data2);
  
  
  mcmc(data1, data2, n1, n2, nb_col1,
       B, dye_swap,
       gamma_1, gamma1_p, gamma_2, gamma2_p,
       mu, mu_p,
       beta2, beta2_p,
       alpha2, alpha2_p,
       delta22,
       delta22_p,
       eta, eta_p,
       lambda_eps1,lambda_eps1_p, lambda_eps2, lambda_eps2_p,
       a_eps,b_eps,
       a_eps_p,b_eps_p,
       w, df_choice, nb_df,
       df, df_p, w_p,
       rho, rho_p, 
       lambda_gamma1, lambda_gamma2,
       lambda_gamma1_p, lambda_gamma2_p,
       min_iter, batch, all_out, verbose);

  PutRNGstate();
  
  free_dmatrix(data1, *n1);
  free_dmatrix(data2, *n1);
  
  
}
 


void mcmc(double **data1,double **data2, int *n1, int *n2, int *nb_col1,
	  int *B, int *dye_swap, double *gamma_1, double *gamma1_p, 
	  double *gamma_2, double *gamma2_p, 
	  double *mu, double *mu_p,
	  double *beta2, double *beta2_p,
	  double *alpha2, double *alpha2_p,
	  double *delta22,
	  double *delta22_p,
	  double *eta, double *eta_p,
	  double *lambda_eps1, double *lambda_eps1_p,
	  double *lambda_eps2, double *lambda_eps2_p,
	  double *a_eps, double *b_eps,
	  double *a_eps_p, double *b_eps_p,
	  double *w, double *df_choice, int *nb_df,
	  double *df, double *df_p, double *w_p,
	  double *rho, double *rho_p,
	  double *lambda_gamma1, double *lambda_gamma2,
	  double *lambda_gamma1_p, double *lambda_gamma2_p,
	  int *min_iter, int *batch, int *all_out, int *verbose)
{
  int i,j,k;
  int count=0,count2=0;

  double lambda_eps1_new, lambda_eps2_new;
  double l_epsilon=0,l_epsilon_new=0;
  double dens_epsilon=0,dens_epsilon_new=0;
  double SSR1=0,SSR2=0,SS12=0;
  double pi=3.14;
  double Sgamma1=0;
  double Sgamma2=0;
  double SSgamma1=0, SSgamma1sd=0, SSgamma2=0, SSgamma2sd=0;
  double SSeps_sd=0.;
  double SSeps_sd_new=0.;
  double *sum_l_w, *sum_l_w_new;
  double dens_t=0,dens_t_new=0;
  double *df_new;
  
  
  /** Parameters for the full conditional of mu **/
  double post_mean_mu=0,post_prec_mu=0;
  double post_mean_beta=0,post_prec_beta=0;
  double post_mean_alpha=0,post_prec_alpha=0;
  double post_mean_delta=0,post_prec_delta=0;
  double post_mean_eta=0,post_prec_eta=0;
  double mean_gamma1=0,mean_gamma2=0;
  

  /** Parameter used in the prior of the gamma's **/
  double *mu1, *mu2;
  

  /** Parameter used to up-date lambda_eps **/
  double sum_res1;
  double sum_res2;

  /** Parameter used in the slice sampling **/
  double sum_lambda=0;
  double sum_log_lambda=0;
  double old_b=0,width_b=5.;
  double old_a=0,width_a=1.;

  /** Parameter used to rescale the gammas **/
  double m_gamma1=0;
  double m_gamma2=0;
  
  /** use for the t-distribution **/
  df_new=dvector(*n2,1);
  sum_l_w=dvector(*n2,0);
  sum_l_w_new=dvector(*n2,0);
  /** Set the prior means to zero **/
  mu1=dvector(1,0);
  mu2=dvector(1,0);
  
  if(*verbose==1)
    printf("--- Model fitting has started --- \n");
  
  for(k=0;k<*B;k++)
    { 
      
      if(*verbose==1 && ((k+1)*100)%(10**B)==0) 
 	{ 
 	  printf("%d percent completed \n",(((k+1)*100)/(*B))); 
	}
      
      
      /** Update the gamma's **/
      Sgamma1=0.;Sgamma2=0.;
      mean_gamma1=0;
      mean_gamma2=0;
      for(i=0;i<*n1;i++)
	{
	  
	  SSgamma1=0;SSgamma2=0;
	  SSgamma1sd=0;SSgamma2sd=0;
	  
	  
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      
	      SSgamma1+=w[j*(*n1)+i]*(lambda_eps1[i]*(data1[i][j]-*mu-eta[j])-*rho*sqrt((lambda_eps1[i]*(lambda_eps2[i])))*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]));
	      SSgamma1sd+=w[j*(*n1)+i]*(lambda_eps1[i]);
	    }
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      
	      SSgamma1+=w[j*(*n1)+i]*(lambda_eps1[i]*(data1[i][j]-*mu-*beta2-eta[j])-*rho*sqrt((lambda_eps1[i]*(lambda_eps2[i])))*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]));
	      SSgamma1sd+=w[j*(*n1)+i]*(lambda_eps1[i]);
	      
	    }
	  
	  gamma_1[i]=rnorm((*lambda_gamma1*(*mu1)+1./(1.-*rho*(*rho))*SSgamma1)/(*lambda_gamma1+1./(1.-*rho*(*rho))*SSgamma1sd),1./sqrt(*lambda_gamma1+1/(1-*rho*(*rho))*SSgamma1sd));
	  /** Use to update b_i in the gibbs **/
	  Sgamma1+=(gamma_1[i]-*mu1)*(gamma_1[i]-*mu1);
	  mean_gamma1+=gamma_1[i];
	  
	  
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      
	      SSgamma2sd+=w[j*(*n1)+i]*(lambda_eps2[i]);
	      SSgamma2+=w[j*(*n1)+i]*(lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-eta[j])-*rho*sqrt((lambda_eps1[i]*(lambda_eps2[i])))*(data1[i][j]-*mu-gamma_1[i]-eta[j]));
	      
	    }
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      
	      SSgamma2sd+=w[j*(*n1)+i]*(lambda_eps2[i]);
	      SSgamma2+=w[j*(*n1)+i]*(lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-eta[j])-*rho*sqrt((lambda_eps1[i]*(lambda_eps2[i])))*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]));
	      
	    }
	  
	  gamma_2[i]=rnorm((*lambda_gamma2*(*mu2)+1./(1.-*rho*(*rho))*SSgamma2)/(*lambda_gamma2+1./(1.-*rho*(*rho))*SSgamma2sd),1/sqrt(*lambda_gamma2+1./(1.-*rho*(*rho))*SSgamma2sd));
	  /** Use to update b_i in the gibbs **/
	  Sgamma2+=(gamma_2[i]-*mu1)*(gamma_2[i]-*mu1);
	  mean_gamma2+=gamma_2[i];
	}
      
      
      /** Up-date lambda_gamma1 and lambda_gamma2 **/
      *lambda_gamma1=rgamma(1.+*n1/2.,1./(0.005+Sgamma1/2.));
      *lambda_gamma2=rgamma(1.+*n1/2.,1./(0.005+Sgamma2/2.));
      
      
            

      /** Up-date the epsilon **/
      for(i=0;i<*n1;i++)
	{
	  

	  SSeps_sd=0.;
	  SSeps_sd_new=0.;
	  
	  SSR1=0;SSR2=0;SS12=0;
	  
	  
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
	      SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      
	  	      
	    }
	  
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
	      SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
	      SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
	      
	  
	    }
	  
	  /** Gibbs sampling if lambda_eps' are independent **/
	  lambda_eps1_new=rgamma(*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR1/2.));

	  lambda_eps2_new=rgamma(*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR2/2.));
	  

	  /** up-date lambda_epsilon **/
	  l_epsilon=0.,l_epsilon_new=0.;
	  SSeps_sd=0.;
	  SSeps_sd_new=0.;
	  
	  SSR1=0;SSR2=0;SS12=0;
	  
	  
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
	      SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      
	      SSeps_sd+=log(w[j*(*n1)+i])+log(lambda_eps1[i])/2.+log(lambda_eps2[i])/2.;
	      SSeps_sd_new+=log(w[j*(*n1)+i])+log(lambda_eps1_new)/2.+log(lambda_eps2_new)/2.;
	      
	    }
	  
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
	      SSR2+=w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
	      SS12+=w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
	      
	      SSeps_sd+=log(w[j*(*n1)+i])+log(lambda_eps1[i])/2.+log(lambda_eps2[i])/2.;
	      SSeps_sd_new+=log(w[j*(*n1)+i])+log(lambda_eps1_new)/2.+log(lambda_eps2_new)/2.;
	      
	    }
	  /** Store another sum (likelihood) to do Metropolis **/
	  l_epsilon+=lambda_eps1[i]*SSR1-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*SS12+lambda_eps2[i]*SSR2;
	  l_epsilon_new+=lambda_eps1_new*SSR1-2*sqrt((lambda_eps1_new)*(lambda_eps2_new))*(*rho)*SS12+lambda_eps2_new*SSR2;
	  
	  
      
	  /** Up-date epsilon1 and epsilon2 by metropolis **/
	  dens_epsilon=-1./2*l_epsilon/(1-*rho*(*rho))+SSeps_sd
	    +dgamma(lambda_eps1[i],*a_eps*(*a_eps)/(*b_eps),(*b_eps)/(*a_eps),1.)+dgamma(lambda_eps2[i],*a_eps*(*a_eps)/(*b_eps),(*b_eps)/(*a_eps),1.);
	  dens_epsilon_new=-1./2*l_epsilon_new/(1-*rho*(*rho))+SSeps_sd_new
	    +dgamma(lambda_eps1_new,*a_eps*(*a_eps)/(*b_eps),(*b_eps)/(*a_eps),1.)+dgamma(lambda_eps2_new,*a_eps*(*a_eps)/(*b_eps),(*b_eps)/(*a_eps),1.);
	  
	  
	  /** Accept the up-dated variables **/
	  if((dens_epsilon_new-dens_epsilon
	      +dgamma(lambda_eps1[i],*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR1/2.),1.)
	      +dgamma(lambda_eps2[i],*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR2/2.),1.)
	      -dgamma(lambda_eps1_new,*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR1/2.),1.)
	      -dgamma(lambda_eps2_new,*a_eps*(*a_eps)/(*b_eps)+*n2/2.,1./(*a_eps/(*b_eps)+SSR2/2.),1.))
	     >log(runif(0,1)))
	    {
	      lambda_eps1[i]=lambda_eps1_new;
	      lambda_eps2[i]=lambda_eps2_new;
	      
	    }
	}
      


      sum_lambda=0;
      sum_log_lambda=0;
      for(i=0;i<*n1;i++)
	{
	  sum_lambda+=lambda_eps1[i]+lambda_eps2[i];
	  sum_log_lambda+=log(lambda_eps1[i])+log(lambda_eps2[i]);
	}

      *a_eps=slice_sampling_a2(*a_eps, width_a, 100, sum_log_lambda, sum_lambda, *b_eps, *n1); 
      *b_eps=slice_sampling_b2(*b_eps, width_b, 100, sum_log_lambda, sum_lambda, *a_eps, *n1); 


      
      SSR1=0;SSR2=0;SS12=0;	  
      /** Up-date rho **/
      for(i=0;i<*n1;i++)
	{
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=lambda_eps1[i]*w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j]);
	      SSR2+=lambda_eps2[i]*w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      SS12+=sqrt(lambda_eps1[i]*lambda_eps2[i])*w[j*(*n1)+i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]);
	      
	    }
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      /** General sum of squares used in the likelihood **/
	      SSR1+=lambda_eps1[i]*w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j]);
	      SSR2+=lambda_eps2[i]*w[j*(*n1)+i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
	      SS12+=sqrt(lambda_eps1[i]*lambda_eps2[i])*w[j*(*n1)+i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]);
     	      
	    }
	  
	}
      *rho=slice_sampling_rho2(*rho, 0.01, 100, SSR1, SSR2, SS12, *n1**n2);
      


      /** Up-date mu **/       
      post_mean_mu=0,post_prec_mu=0;

      for(i=0;i<*n1;i++)
	{
	  for(j=0;j<(*nb_col1);j++)
	    {
	      post_prec_mu+=w[j*(*n1)+i]*(lambda_eps1[i]+lambda_eps2[i]-2*(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i])))/(1.-(*rho)*(*rho));
	      post_mean_mu+=w[j*(*n1)+i]*(-(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i]))*(data1[i][j]-gamma_1[i]-eta[j]+data2[i][j]-*alpha2-gamma_2[i]-eta[j])+(lambda_eps1[i]*(data1[i][j]-gamma_1[i]-eta[j])+lambda_eps2[i]*(data2[i][j]-*alpha2-gamma_2[i]-eta[j])));
	      
	    }
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_prec_mu+=w[j*(*n1)+i]*(lambda_eps1[i]+lambda_eps2[i]-2*(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i])))/(1.-(*rho)*(*rho));
	      post_mean_mu+=w[j*(*n1)+i]*(-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*beta2-gamma_1[i]-eta[j]+data2[i][j]-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])+(lambda_eps1[i])*(data1[i][j]-*beta2-gamma_1[i]-eta[j])+(lambda_eps2[i])*(data2[i][j]-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]));
	    }
	}
      post_prec_mu+=0.04;
      post_mean_mu=(1./(1.-(*rho)*(*rho)))*post_mean_mu/post_prec_mu;
      *mu=rnorm(post_mean_mu,1./sqrt(post_prec_mu));
	 
      /** If dye swap experiment **/
      if(*dye_swap==1 && (*nb_col1>0))
	{
	  /** Up-date beta **/
	  post_mean_beta=0,post_prec_beta=0;
	  
	  for(i=0;i<*n1;i++)
	    {
	      for(j=(*nb_col1);j<*n2;j++)
		{
		  post_prec_beta+=w[j*(*n1)+i]*(lambda_eps1[i]+lambda_eps2[i]-2*(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i])))/(1.-(*rho)*(*rho));
		  post_mean_beta+=w[j*(*n1)+i]*(-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*mu-gamma_1[i]-eta[j]+data2[i][j]-*mu-*alpha2-*delta22-gamma_2[i]-eta[j])+(lambda_eps1[i])*(data1[i][j]-*mu-gamma_1[i]-eta[j])+(lambda_eps2[i])*(data2[i][j]-*mu-*alpha2-*delta22-gamma_2[i]-eta[j]));
		}
	    }
	  
	  post_prec_beta+=0.04;
	  post_mean_beta=(1./(1.-(*rho)*(*rho)))*post_mean_beta/post_prec_beta;
	  *beta2=rnorm(post_mean_beta,1./sqrt(post_prec_beta));
	  
	  
	  /** Up-date delta **/
	  post_mean_delta=0,post_prec_delta=0;
	  
	  for(i=0;i<*n1;i++)
	    {
	      
	      for(j=(*nb_col1);j<*n2;j++)
		{
		  
		  post_prec_delta+=w[j*(*n1)+i]*(lambda_eps2[i])/(1.-(*rho)*(*rho));
		  post_mean_delta+=w[j*(*n1)+i]*(-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])+(lambda_eps2[i])*(data2[i][j]-*mu-*alpha2-*beta2-gamma_2[i]-eta[j]));
		}
	    }
	  
	  post_prec_delta+=0.04;
	  post_mean_delta=(1./(1.-(*rho)*(*rho)))*post_mean_delta/post_prec_delta;
	  *delta22=rnorm(post_mean_delta,1./sqrt(post_prec_delta));
	  
	}
      
      
      /** Up-date eta **/       
      if((*nb_col1>=1 && (*dye_swap)==0) | (*nb_col1>=2 && (*dye_swap)==1))
	{
	  eta[0]=0;
	  for(j=1;j<(*nb_col1);j++)
	    {
	      post_mean_eta=0,post_prec_eta=0;
	      for(i=0;i<*n1;i++)
		{
		  post_prec_eta+=w[j*(*n1)+i]*(lambda_eps1[i]+lambda_eps2[i]-2*(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i])))/(1.-(*rho)*(*rho));
		  post_mean_eta+=w[j*(*n1)+i]*((lambda_eps1[i])*(data1[i][j]-*mu-gamma_1[i])+(lambda_eps2[i])*(data2[i][j]-*mu-*alpha2-gamma_2[i])-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*mu-gamma_1[i]+data2[i][j]-*mu-*alpha2-gamma_2[i]));
		}
	      post_prec_eta+=0.04;
	      post_mean_eta=(1./(1.-(*rho)*(*rho)))*post_mean_eta/post_prec_eta;
	      eta[j]=rnorm(post_mean_eta,1./sqrt(post_prec_eta));
	      
	    }
	  
	  
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_mean_eta=0,post_prec_eta=0;
	      for(i=0;i<*n1;i++)
		{
		  post_prec_eta+=w[j*(*n1)+i]*(lambda_eps1[i]+lambda_eps2[i]-2*(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i])))/(1.-(*rho)*(*rho));
		  post_mean_eta+=w[j*(*n1)+i]*((lambda_eps1[i])*(data1[i][j]-*mu-*beta2-gamma_1[i])+(lambda_eps2[i])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i])-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*mu-*beta2-gamma_1[i]+data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]));
		}
	      post_prec_eta+=0.04;
	      post_mean_eta=(1./(1.-(*rho)*(*rho)))*post_mean_eta/post_prec_eta;
	      eta[j]=rnorm(post_mean_eta,1./sqrt(post_prec_eta));	      
	    }
	  if(*dye_swap==0)
	    eta[*n2-1]=0.;
	}
      else
	{
	  /** Not enough data to estimate the array effect **/
	  for(j=0;j<*n2;j++)
	    eta[j]=0;
	}
      
      

      /** Up-date alpha **/
      post_mean_alpha=0,post_prec_alpha=0;

      for(i=0;i<*n1;i++)
	{
	  for(j=0;j<(*nb_col1);j++)
	    {
	      post_prec_alpha+=w[j*(*n1)+i]*(lambda_eps2[i])/(1.-(*rho)*(*rho));
	      post_mean_alpha+=w[j*(*n1)+i]*(-(*rho)*sqrt(lambda_eps1[i]*(lambda_eps2[i]))*(data1[i][j]-*mu-gamma_1[i]-eta[j])+(lambda_eps2[i]*(data2[i][j]-*mu-gamma_2[i]-eta[j])));
	    }
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_prec_alpha+=w[j*(*n1)+i]*(lambda_eps2[i])/(1.-(*rho)*(*rho));
	      post_mean_alpha+=w[j*(*n1)+i]*(-(*rho)*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])+(lambda_eps2[i])*(data2[i][j]-*mu-*delta22-*beta2-gamma_2[i]-eta[j]));
	    }
	}
      
      post_prec_alpha+=0.04;
      post_mean_alpha=(1./(1.-(*rho)*(*rho)))*post_mean_alpha/post_prec_alpha;
      *alpha2=rnorm(post_mean_alpha,1./sqrt(post_prec_alpha));
      
      
      /** Up-date the weight and the df for the t-distribution **/
      /** A block up-date is used **/
      /** First the parameter w is integrated out and df is up-dated **/
      /** Then the w are up-dated by Gibbs sampling **/

	for(j=0;j<*n2;j++)
	  {
	    
	    /** Up-date the degrees of freedom **/
	    df_new[j]=df_choice[(int)(runif(0.,1.)*(*nb_df))];  
	    sum_l_w[j]=0.;
	    sum_l_w_new[j]=0.;
	  }
      

	
	for(i=0;i<(*n1);i++)
	  { 
	    for(j=0;j<(*nb_col1);j++)
	      {
		
		sum_l_w[j]+=log(2./df[j])-(df[j]/2.+1)*log(1+1/((1-*rho*(*rho)))
							   *(lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
							     -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)
							     *(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
							     +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))/df[j])
		  +lgammafn(df[j]/2.+1)-lgammafn(df[j]/2.);
		
		sum_l_w_new[j]+=log(2./df_new[j])-(df_new[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
								       (lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
									-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-gamma_1[i]-eta[j])
									*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
									+lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))/df_new[j])
		  +lgammafn(df_new[j]/2.+1)-lgammafn(df_new[j]/2.);
		
	      }
	  
	    for(j=(*nb_col1);j<*n2;j++)
	      {
	      
		sum_l_w[j]+=log(2./df[j])-(df[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
							   (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
							    -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])
							    +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]))/df[j])+lgammafn(df[j]/2.+1)-lgammafn(df[j]/2.);
		
		
		sum_l_w_new[j]+=log(2./df_new[j])-(df_new[j]/2.+1)*log(1+1/((1-*rho*(*rho)))*
								       (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
									-2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])
									+lambda_eps2[i]*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]))/df_new[j])+lgammafn(df_new[j]/2.+1)-lgammafn(df_new[j]/2.);
	      }
	  }
	
	
	/** Up-date the degrees of freedom for the t-distribution **/
	
	for(j=0;j<(*nb_col1);j++)
	  {
	    if((sum_l_w_new[j]-sum_l_w[j])>log(runif(0,1)))
	      {
		df[j]=df_new[j];
	      }
	    for(i=0;i<*n1;i++)
	      w[j*(*n1)+i]=rgamma(df[j]/2.+1,1./(df[j]/2.+1/(2*(1-*rho*(*rho)))*
						 (lambda_eps1[i]*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data1[i][j]-*mu-gamma_1[i]-eta[j])
						  -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])
						  +lambda_eps2[i]*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*alpha2-gamma_2[i]-eta[j]))));
	  }
	
	for(j=(*nb_col1);j<*n2;j++)
	  {
	    if((sum_l_w_new[j]-sum_l_w[j])>log(runif(0,1)))
	      {
		df[j]=df_new[j];
	      }
	    for(i=0;i<*n1;i++)
	      w[j*(*n1)+i]=rgamma(df[j]/2.+1,1./(df[j]/2.+1/(2*(1-*rho*(*rho)))*
						 (lambda_eps1[i]*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])
						  -2*sqrt((lambda_eps1[i])*(lambda_eps2[i]))*(*rho)*(data1[i][j]-*mu-*beta2-gamma_1[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])
						  +lambda_eps2[i]*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])*(data2[i][j]-*mu-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]))));
	  }
	
	


	
	
	if(k>=*min_iter)
	  {
	    if(((count2+1)%(*batch))==0) /** Batch sampling **/
	      {
		
		/** Mean of the gamma parameters **/
		m_gamma1=mean_vec(gamma_1, n1);
		m_gamma2=mean_vec(gamma_2, n1);
		  
		for(i=0;i<*n1;i++)
		  {
		    
		    if(*all_out==1)
		      {
			for(j=0;j<*n2;j++)
			  {
			    /** Compute the mean of the weights **/
			    w_p[j*(*n1)+i]+=w[j*(*n1)+i]/((*B-*min_iter)/(*batch));		  
			    df_p[count*(*n2)+j]=df[j];
			    eta_p[count*(*n2)+j]=eta[j];
			  }
			lambda_eps1_p[count*(*n1)+i]=lambda_eps1[i];
			lambda_eps2_p[count*(*n1)+i]=lambda_eps2[i];
		      
			gamma1_p[count*(*n1)+i]=gamma_1[i]-m_gamma1;
			gamma2_p[count*(*n1)+i]=gamma_2[i]-m_gamma2;
		      }
		    else /** Posterior mean **/
		      {
		      	for(j=0;j<*n2;j++)
			  {
			    /** Compute the mean of the weights **/
			    w_p[j*(*n1)+i]+=w[j*(*n1)+i]/((*B-*min_iter)/(*batch));		  
			    df_p[j]+=df[j]/((*B-*min_iter)/(*batch));
			    eta_p[j]+=eta[j]/((*B-*min_iter)/(*batch));
			  }
			lambda_eps1_p[i]+=lambda_eps1[i]/((*B-*min_iter)/(*batch));
			lambda_eps2_p[i]+=lambda_eps2[i]/((*B-*min_iter)/(*batch));
			
			gamma1_p[i]+=(gamma_1[i]-m_gamma1)/((*B-*min_iter)/(*batch));
			gamma2_p[i]+=(gamma_2[i]-m_gamma2)/((*B-*min_iter)/(*batch));
		      }
		  }
		

		
		if(*all_out==1)
		  {
		    mu_p[count]=*mu+m_gamma1;
		    beta2_p[count]=*beta2;
		    alpha2_p[count]=*alpha2+m_gamma2-m_gamma1;
		    delta22_p[count]=*delta22;
		    a_eps_p[count]=*a_eps;
		    b_eps_p[count]=*b_eps;
		    lambda_gamma1_p[count]=*lambda_gamma1;
		    lambda_gamma2_p[count]=*lambda_gamma2;
		    rho_p[count]=*rho;
		  }
		else /** Compute the posterior mean **/
		  {
		    *mu_p+=(*mu+m_gamma1)/((*B-*min_iter)/(*batch));
		    *beta2_p+=*beta2/((*B-*min_iter)/(*batch));
		    *alpha2_p+=(*alpha2+m_gamma2-m_gamma1)/((*B-*min_iter)/(*batch));
		    *delta22_p+=*delta22/((*B-*min_iter)/(*batch));
		    *a_eps_p+=*a_eps/((*B-*min_iter)/(*batch));
		    *b_eps_p+=*b_eps/((*B-*min_iter)/(*batch));
		    *lambda_gamma1_p+=*lambda_gamma1/((*B-*min_iter)/(*batch));
		    *lambda_gamma2_p+=*lambda_gamma2/((*B-*min_iter)/(*batch));
		    *rho_p+=*rho/((*B-*min_iter)/(*batch));
		  }
		
		count++;
	      }
	    count2++;
	  }
	
	
	
    }

  
  Free(df_new);
  Free(sum_l_w);
  Free(sum_l_w_new);
  Free(mu1);
  Free(mu2);
  
  
} 
