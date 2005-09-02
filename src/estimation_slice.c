#include "util.h" 

void mcmc_shift(double **data1,double **data2, int *n1, int *n2, int *nb_col1,
		int *dye_swap,
		int *B, double *gamma_1, double *gamma1_p, 
		double *gamma_2, double *gamma2_p, 
		double *mu, double *mu_p,
		double *beta2, double *beta2_p,
		double *alpha2, double *alpha2_p,
		double *delta22,
		double *delta22_p,
		double *eta, double *eta_p,
		double *lambda_eps1, double *lambda_eps1_p,
		double *lambda_eps2, double *lambda_eps2_p,
		double *w, double *df_choice, int *nb_df,
		double *df, 
		double *rho, double *rho_p,
		double *shift, double *shift_p, 
		double *min_shift,
		double *lambda_gamma1, double *lambda_gamma2,
		double *lambda_gamma1_p, double *lambda_gamma2_p,
		int *min_iter, int *batch, int *all_out);

double log_f_lambda_eps(double SSR1, double SSR2, double SS12, double rho, int n, double lambda_eps1, double lambda_eps2, double a_eps, double b_eps);
void up_date_error_precisions_slice(double **data1, double **data2, int n1, int n2, int nb_col1, double shift, double mu, double alpha2, double beta2, double delta22, double *eta, double *gamma_1, double *gamma_2, double rho, double *lambda_eps1, double *lambda_eps2, double a_eps, double b_eps, double *w);
double slice_sampling_lambda_eps(double width, int p, double SSR1, double SSR2, double SS12, double rho, int n, double lambda_eps1, double lambda_eps2, double a_eps, double b_eps);

void R_link_mcmc_shift(double *data_vec1, double *data_vec2, 
		       int *n1, int *n2, int *nb_col1, int *dye_swap,
		       int *B, 
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
		       double *w, double *df_choice, int *nb_df,
		       double *df, 
		       double *rho, double *rho_p,
		       double *shift, double *shift_p, 
		       double *min_shift,
		       double *lambda_gamma1, double *lambda_gamma2,
		       double *lambda_gamma1_p, double *lambda_gamma2_p,
		       int *min_iter, int *batch, int *all_out)
{
  
  double **data1;
  double **data2;
  
  GetRNGstate();
  data1=dmatrix(*n1, *n2);  
  vec_mat(data_vec1,n1,n2,data1);
  data2=dmatrix(*n1, *n2);  
  vec_mat(data_vec2,n1,n2,data2);

  
  mcmc_shift(data1, data2, n1, n2, nb_col1, dye_swap,
	     B, gamma_1, gamma1_p, gamma_2, gamma2_p,
	     mu, mu_p,
	     beta2, beta2_p,
	     alpha2, alpha2_p,
	     delta22,
	     delta22_p,
	     eta, eta_p,
	     lambda_eps1,lambda_eps1_p, 
	     lambda_eps2, lambda_eps2_p,
	     w, df_choice, nb_df,
	     df,
	     rho, rho_p,
	     shift, shift_p,
	     min_shift,
	     lambda_gamma1, lambda_gamma2,
	     lambda_gamma1_p, lambda_gamma2_p,
	     min_iter, batch, all_out);

  PutRNGstate();
  
  free_dmatrix(data1, *n1);
  free_dmatrix(data2, *n1);
  
  
}
 


void mcmc_shift(double **data1,double **data2, int *n1, int *n2, int *nb_col1, int *dye_swap,
		int *B, double *gamma_1, double *gamma1_p, 
		double *gamma_2, double *gamma2_p, 
		double *mu, double *mu_p,
		double *beta2, double *beta2_p,
		double *alpha2, double *alpha2_p,
		double *delta22,
		double *delta22_p,
		double *eta, double *eta_p,
		double *lambda_eps1, double *lambda_eps1_p,
		double *lambda_eps2, double *lambda_eps2_p,
		double *w, double *df_choice, int *nb_df,
		double *df,
		double *rho, double *rho_p,
		double *shift, double *shift_p, 
		double *min_shift,
		double *lambda_gamma1, double *lambda_gamma2,
		double *lambda_gamma1_p, double *lambda_gamma2_p,
		int *min_iter, int *batch, int *all_out)
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
  double *df_new;
  
  
  
  /** Parameters for the full conditional of mu **/
  double post_mean_mu=0,post_prec_mu=0;
  double post_mean_beta=0,post_prec_beta=0;
  double post_mean_alpha=0,post_prec_alpha=0;
  double post_mean_delta=0,post_prec_delta=0;
  double post_mean_eta=0,post_prec_eta=0;
  double mean_gamma1=0,mean_gamma2=0;
  double *mu1,*mu2;
  double m_gamma1=0, m_gamma2=0;


  /** use for the t-distribution **/
  df_new=dvector(*n2,1);
  sum_l_w=dvector(*n2,0);
  sum_l_w_new=dvector(*n2,0);
  mu1=dvector(1,0);
  mu2=dvector(1,0);
  
  
  for(k=0;k<*B;k++)
    {
      
/*       if(((k+1)*100)%(10**B)==0) */
/* 	{ */
/* 	  printf("%d percent completed \n",(((k+1)*100)/(*B))); */
/* 	} */
      
      

      *shift=slice_sampling_shift2(*shift, 2, 200, data1, data2, n1, n2, nb_col1,
				   gamma_1, gamma_2, 
				   mu, beta2, 
				   alpha2, delta22,
				   eta, lambda_eps1,  
				   lambda_eps2, w, rho);
      
      
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
	      
	      SSgamma1+=w[j*(*n1)+i]*(*lambda_eps1*(log2(data1[i][j]+*shift)-*mu-eta[j])-*rho*sqrt((*lambda_eps1*(*lambda_eps2)))*(log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i]-eta[j]));
	      SSgamma1sd+=w[j*(*n1)+i]*(*lambda_eps1);
	    }
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      
	      SSgamma1+=w[j*(*n1)+i]*(*lambda_eps1*(log2(data1[i][j]+*shift)-*mu-*beta2-eta[j])-*rho*sqrt((*lambda_eps1*(*lambda_eps2)))*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]));
	      SSgamma1sd+=w[j*(*n1)+i]*(*lambda_eps1);
	      
	    }
	  
	  gamma_1[i]=rnorm((*lambda_gamma1*(*mu1)+1./(1.-*rho*(*rho))*SSgamma1)/(*lambda_gamma1+1./(1.-*rho*(*rho))*SSgamma1sd),1./sqrt(*lambda_gamma1+1/(1-*rho*(*rho))*SSgamma1sd));

	  /** Use to update b_i in the gibbs **/
	  Sgamma1+=(gamma_1[i]-*mu1)*(gamma_1[i]-*mu1);
	  mean_gamma1+=gamma_1[i];
	  
	  
	  /** First part of the loop before dye swap **/
	  for(j=0;j<(*nb_col1);j++)
	    {
	      
	      SSgamma2sd+=w[j*(*n1)+i]*(*lambda_eps2);
	      SSgamma2+=w[j*(*n1)+i]*(*lambda_eps2*(log2(data2[i][j]+*shift)-*mu-*alpha2-eta[j])-*rho*sqrt((*lambda_eps1*(*lambda_eps2)))*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j]));
	      
	    }
	  
	  /** Second part of the loop after dye swap **/
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      
	      SSgamma2sd+=w[j*(*n1)+i]*(*lambda_eps2);
	      SSgamma2+=w[j*(*n1)+i]*(*lambda_eps2*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-eta[j])-*rho*sqrt((*lambda_eps1*(*lambda_eps2)))*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j]));
	      
	    }
	  
	  gamma_2[i]=rnorm((*lambda_gamma2*(*mu2)+1./(1.-*rho*(*rho))*SSgamma2)/(*lambda_gamma2+1./(1.-*rho*(*rho))*SSgamma2sd),1/sqrt(*lambda_gamma2+1./(1.-*rho*(*rho))*SSgamma2sd));
	  /** Use to update b_i in the gibbs **/
	  Sgamma2+=(gamma_2[i]-*mu1)*(gamma_2[i]-*mu1);
	  mean_gamma2+=gamma_2[i];
	}
      

      
      /** Up-date lambda_gamma1 and lambda_gamma2 **/
      *lambda_gamma1=rgamma(1.+*n1/2.,1./(0.005+Sgamma1/2.));
      *lambda_gamma2=rgamma(1.+*n1/2.,1./(0.005+Sgamma2/2.));
      
      
      up_date_error_precisions_slice(data1, data2, *n1, *n2, *nb_col1, *shift, *mu, *alpha2, *beta2, *delta22, eta, gamma_1, gamma_2, *rho, lambda_eps1, lambda_eps2, 200., 40000., w);
      
      
      

/*       SSR1=0;SSR2=0;SS12=0; */
/*       for(i=0;i<*n1;i++) */
/* 	{ */
	  
	  

/* 	  for(j=0;j<(*nb_col1);j++) */
/* 	    { */

/* 	      SSR1+=*lambda_eps1*w[j*(*n1)+i]*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j])*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j]); */
/* 	      SSR2+=*lambda_eps2*w[j*(*n1)+i]*(log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i]-eta[j])*(log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i]-eta[j]); */
/* 	      SS12+=sqrt(*lambda_eps1**lambda_eps2)*w[j*(*n1)+i]*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j])*(log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i]-eta[j]); */
	      
/* 	    } */
	  

/* 	  for(j=(*nb_col1);j<*n2;j++) */
/* 	    { */

/* 	      SSR1+=*lambda_eps1*w[j*(*n1)+i]*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j])*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j]); */
/* 	      SSR2+=*lambda_eps1*w[j*(*n1)+i]*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j])*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]); */
/* 	      SS12+=sqrt(*lambda_eps1**lambda_eps2)*w[j*(*n1)+i]*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j])*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-gamma_2[i]-eta[j]); */
     	      
/* 	    } */
/* 	} */
      *rho=0;
      //*rho=slice_sampling_rho2(*rho, 0.01, 10, SSR1, SSR2, SS12, *n1**n2);


      /** Up-date mu **/       
      post_mean_mu=0,post_prec_mu=0;

      for(i=0;i<*n1;i++)
	{
	  for(j=0;j<(*nb_col1);j++)
	    {
	      post_prec_mu+=w[j*(*n1)+i];
	      post_mean_mu+=w[j*(*n1)+i]*(-(*rho)*sqrt(*lambda_eps1*(*lambda_eps2))*(log2(data1[i][j]+*shift)-gamma_1[i]-eta[j]+log2(data2[i][j]+*shift)-*alpha2-gamma_2[i]-eta[j])+(*lambda_eps1*(log2(data1[i][j]+*shift)-gamma_1[i]-eta[j])+*lambda_eps2*(log2(data2[i][j]+*shift)-*alpha2-gamma_2[i]-eta[j])));
	      
	    }
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_prec_mu+=w[j*(*n1)+i];
	      post_mean_mu+=w[j*(*n1)+i]*(-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*beta2-gamma_1[i]-eta[j]+log2(data2[i][j]+*shift)-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j])+(*lambda_eps1)*(log2(data1[i][j]+*shift)-*beta2-gamma_1[i]-eta[j])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*beta2-*alpha2-*delta22-gamma_2[i]-eta[j]));
	    }
	}
      post_prec_mu=post_prec_mu*(*lambda_eps1+*lambda_eps2-2*(*rho)*sqrt(*lambda_eps1*(*lambda_eps2)))/(1.-(*rho)*(*rho))+0.04;
      post_mean_mu=(1./(1.-(*rho)*(*rho)))*post_mean_mu/post_prec_mu;
      *mu=rnorm(post_mean_mu,1./sqrt(post_prec_mu));
	  
      
      if(*dye_swap==1 && (*nb_col1>0))
	{
	  /** Up-date beta **/
	  post_mean_beta=0,post_prec_beta=0;
	  
	  for(i=0;i<*n1;i++)
	    { 
	      for(j=(*nb_col1);j<*n2;j++)
		{
		  post_prec_beta+=w[j*(*n1)+i];
		  post_mean_beta+=w[j*(*n1)+i]*(-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j]+log2(data2[i][j]+*shift)-*mu-*alpha2-*delta22-gamma_2[i]-eta[j])+(*lambda_eps1)*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*mu-*alpha2-*delta22-gamma_2[i]-eta[j]));
		}
	    }
	  
	  post_prec_beta=post_prec_beta*(*lambda_eps1+*lambda_eps2-2*(*rho)*sqrt(*lambda_eps1*(*lambda_eps2)))/(1.-(*rho)*(*rho))+0.04;
	  post_mean_beta=(1./(1.-(*rho)*(*rho)))*post_mean_beta/post_prec_beta;
	  *beta2=rnorm(post_mean_beta,1./sqrt(post_prec_beta));
	  
	  /** Up-date delta **/
	  post_mean_delta=0,post_prec_delta=0;
	  for(i=0;i<*n1;i++)
	    {
	      
	      for(j=(*nb_col1);j<*n2;j++)
		{
		  
		  post_prec_delta+=w[j*(*n1)+i];
		  post_mean_delta+=w[j*(*n1)+i]*(-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-gamma_2[i]-eta[j]));
		}
	    }
      
	  post_prec_delta=post_prec_delta*(*lambda_eps2)/(1.-(*rho)*(*rho))+0.04;
	  post_mean_delta=(1./(1.-(*rho)*(*rho)))*post_mean_delta/post_prec_delta;
	  *delta22=rnorm(post_mean_delta,1./sqrt(post_prec_delta));
	}


      /** Up-date eta **/    
      if((*nb_col1>=1 && (*dye_swap)==0) | (*nb_col1>=2 && (*dye_swap)==1))
	{
	  eta[0]=0.;
	  for(j=1;j<(*nb_col1);j++)
	    {
	      post_mean_eta=0,post_prec_eta=0;
	      for(i=0;i<*n1;i++)
		{
		  post_prec_eta+=w[j*(*n1)+i];
		  post_mean_eta+=w[j*(*n1)+i]*((*lambda_eps1)*(log2(data1[i][j]+*shift)-*mu-gamma_1[i])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i])-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]+log2(data2[i][j]+*shift)-*mu-*alpha2-gamma_2[i]));
		}
	      post_prec_eta=post_prec_eta*(*lambda_eps1+*lambda_eps2-2*(*rho)*sqrt(*lambda_eps1*(*lambda_eps2)))/(1.-(*rho)*(*rho))+0.04;
	      post_mean_eta=(1./(1.-(*rho)*(*rho)))*post_mean_eta/post_prec_eta;
	      eta[j]=rnorm(post_mean_eta,1./sqrt(post_prec_eta));
	      
	    }
	  
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_mean_eta=0,post_prec_eta=0;
	      for(i=0;i<*n1;i++)
		{
		  post_prec_eta+=w[j*(*n1)+i];
		  post_mean_eta+=w[j*(*n1)+i]*((*lambda_eps1)*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*mu-*beta2-*alpha2-*delta22-gamma_2[i])-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]+log2(data2[i][j]+*shift)-*mu-*alpha2-*beta2-*delta22-gamma_2[i]));
		}
	      post_prec_eta=post_prec_eta*(*lambda_eps1+*lambda_eps2-2*(*rho)*sqrt(*lambda_eps1*(*lambda_eps2)))/(1.-(*rho)*(*rho))+0.04;
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
	      post_prec_alpha+=w[j*(*n1)+i];
	      post_mean_alpha+=w[j*(*n1)+i]*(-(*rho)*sqrt(*lambda_eps1*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-gamma_1[i]-eta[j])+(*lambda_eps2*(log2(data2[i][j]+*shift)-*mu-gamma_2[i]-eta[j])));
	    }
	  for(j=(*nb_col1);j<*n2;j++)
	    {
	      post_prec_alpha+=w[j*(*n1)+i];
	      post_mean_alpha+=w[j*(*n1)+i]*(-(*rho)*sqrt((*lambda_eps1)*(*lambda_eps2))*(log2(data1[i][j]+*shift)-*mu-*beta2-gamma_1[i]-eta[j])+(*lambda_eps2)*(log2(data2[i][j]+*shift)-*mu-*delta22-*beta2-gamma_2[i]-eta[j]));
	    }
	}
      
      post_prec_alpha=post_prec_alpha*(*lambda_eps2)/(1.-(*rho)*(*rho))+0.04;
      post_mean_alpha=(1./(1.-(*rho)*(*rho)))*post_mean_alpha/post_prec_alpha;
      *alpha2=rnorm(post_mean_alpha,1./sqrt(post_prec_alpha));
      

	if(k>=*min_iter)
	  {
	    if(((count2+1)%(*batch))==0) /** Batch sampling **/
	      {
		/** Mean of the gamma parameters **/
		m_gamma1=mean_vec(gamma_1, n1);
		m_gamma2=mean_vec(gamma_2, n1);
		if(*all_out==1)
		  {
		    for(i=0;i<*n1;i++)
		      {
			for(j=0;j<*n2;j++)
			  {
			    eta_p[count*(*n2)+j]=eta[j];
			  }
			gamma1_p[count*(*n1)+i]=gamma_1[i]-m_gamma1;
			gamma2_p[count*(*n1)+i]=gamma_2[i]-m_gamma2;
		      }
		    shift_p[count]=*shift;
		    mu_p[count]=*mu+m_gamma1;
		    beta2_p[count]=*beta2;
		    alpha2_p[count]=*alpha2+m_gamma2-m_gamma1;
		    delta22_p[count]=*delta22;
		
		    lambda_eps1_p[count]=*lambda_eps1;
		    lambda_eps2_p[count]=*lambda_eps2;
		    lambda_gamma1_p[count]=*lambda_gamma1;
		    lambda_gamma2_p[count]=*lambda_gamma2;
		    rho_p[count]=*rho;
		  }
		else /** Posterior mean **/
		  {
		    for(i=0;i<*n1;i++)
		      {
			for(j=0;j<*n2;j++)
			  {
			    eta_p[j]+=eta[j]/((*B-*min_iter)/(*batch));
			  }
			gamma1_p[i]+=(gamma_1[i]-m_gamma1)/((*B-*min_iter)/(*batch));
			gamma2_p[i]+=(gamma_2[i]-m_gamma2)/((*B-*min_iter)/(*batch));
		      }
		    *shift_p+=*shift/((*B-*min_iter)/(*batch));
		    *mu_p+=(*mu+m_gamma1)/((*B-*min_iter)/(*batch));
		    *beta2_p+=*beta2/((*B-*min_iter)/(*batch));
		    *alpha2_p+=(*alpha2+m_gamma2-m_gamma1)/((*B-*min_iter)/(*batch));
		    *delta22_p+=*delta22/((*B-*min_iter)/(*batch));
		
		    *lambda_eps1_p+=*lambda_eps1/((*B-*min_iter)/(*batch));
		    *lambda_eps2_p+=*lambda_eps2/((*B-*min_iter)/(*batch));
		    *lambda_gamma1_p+=*lambda_gamma1/((*B-*min_iter)/(*batch));
		    *lambda_gamma2_p+=*lambda_gamma2/((*B-*min_iter)/(*batch));
		    *rho_p+=*rho/((*B-*min_iter)/(*batch));
		  }
		
		count++;
	      }
	    count2++;
	  }
	
	
	
    }

  free(mu1);
  free(mu2);
  free(df_new);
  free(sum_l_w);
  free(sum_l_w_new);
  
   
}

double slice_sampling_lambda_eps(double width, int p, double SSR1, double SSR2, double SS12, double rho, int n, double lambda_eps1, double lambda_eps2, double a_eps, double b_eps)
{

  double L, R, z, lambda_eps_new;
  int J,K;
  double log_f_R=0,log_f_L=0;


  z=log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, lambda_eps1, lambda_eps2,a_eps,b_eps)-rgamma(1,1);

  /** lower bound **/
  L=lambda_eps1-width*runif(0,1);
  R=L+width;
  J=(int)(p*runif(0,1));
  K=(p-1)-J;
  

  /** Find the interval **/
  log_f_L=log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, L, lambda_eps2,a_eps,b_eps);
  log_f_R=log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, R, lambda_eps2,a_eps,b_eps);

  /** Find the interval **/
  while(J>0 & z<log_f_L)
    {
      
      L=L-width;
      log_f_L=log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, L, lambda_eps2,a_eps,b_eps);
      
      J--;
    }
  
  while(K>0 & z<log_f_R)
    {
      
      R=R+width;
      log_f_R=log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, R, lambda_eps2,a_eps,b_eps);
      
      K--;
    }
  
  
  /** Make sure it is in the range **/  
  L=fmax2(0.,L);
  
  lambda_eps_new=runif(L,R);
  while(log_f_lambda_eps(SSR1, SSR2, SS12, rho, n, lambda_eps_new, lambda_eps2,a_eps,b_eps)<z)
    {
      if(lambda_eps_new<lambda_eps1)
      	L=lambda_eps_new;
      else
      	R=lambda_eps_new;
      lambda_eps_new=runif(L,R);
    }
  
  return(lambda_eps_new);
  
}

double log_f_lambda_eps(double SSR1, double SSR2, double SS12, double rho, int n, double lambda_eps1, double lambda_eps2, double a_eps, double b_eps)
{
  
  double log_f;
  log_f=n/2.*log(lambda_eps1*lambda_eps2)-1./(2.*(1.-rho*rho))*(lambda_eps1*SSR1-2*rho*sqrt(lambda_eps1*lambda_eps2)*SS12+lambda_eps2*SSR2)
    +dgamma(lambda_eps1,a_eps*a_eps/b_eps,b_eps/a_eps,1.)+dgamma(lambda_eps2,a_eps*a_eps/b_eps,b_eps/a_eps,1.);
  return(log_f);
  
}

void up_date_error_precisions_slice(double **data1, double **data2, int n1, int n2, int nb_col1, double shift, double mu, double alpha2, double beta2, double delta22, double *eta, double *gamma_1, double *gamma_2, double rho, double *lambda_eps1, double *lambda_eps2, double a_eps, double b_eps, double *w)
{
  int i,j;
  double SSR1,SSR2,SS12;

  
  SSR1=0;SSR2=0;SS12=0;
  for(i=0;i<n1;i++)
    {
      /** First part of the loop before dye swap **/
      for(j=0;j<(nb_col1);j++)
	{
	  /** General sum of squares used in the likelihood **/
	  SSR1+=w[j*n1+i]*(log2(data1[i][j]+shift)-mu-gamma_1[i]-eta[j])*(log2(data1[i][j]+shift)-mu-gamma_1[i]-eta[j]);
	  SSR2+=w[j*n1+i]*(log2(data2[i][j]+shift)-mu-alpha2-gamma_2[i]-eta[j])*(log2(data2[i][j]+shift)-mu-alpha2-gamma_2[i]-eta[j]);
	  SS12+=w[j*n1+i]*(log2(data1[i][j]+shift)-mu-gamma_1[i]-eta[j])*(log2(data2[i][j]+shift)-mu-alpha2-gamma_2[i]-eta[j]);
	      
	}
	  
      /** Second part of the loop after dye swap **/
      for(j=nb_col1;j<n2;j++)
	{
	  /** General sum of squares used in the likelihood **/
	  SSR1+=w[j*n1+i]*(log2(data1[i][j]+shift)-mu-beta2-gamma_1[i]-eta[j])*(log2(data1[i][j]+shift)-mu-beta2-gamma_1[i]-eta[j]);
	  SSR2+=w[j*n1+i]*(log2(data2[i][j]+shift)-mu-alpha2-beta2-delta22-gamma_2[i]-eta[j])*(log2(data2[i][j]+shift)-mu-alpha2-beta2-delta22-gamma_2[i]-eta[j]);
	  SS12+=w[j*n1+i]*(log2(data1[i][j]+shift)-mu-beta2-gamma_1[i]-eta[j])*(log2(data2[i][j]+shift)-mu-alpha2-beta2-delta22-gamma_2[i]-eta[j]);
	      
	}
    }
  
  *lambda_eps1=slice_sampling_lambda_eps(.1, 10, SSR1, SSR2, SS12, rho, n1*n2, *lambda_eps1, *lambda_eps2, a_eps, b_eps);
  *lambda_eps2=slice_sampling_lambda_eps(.1, 10, SSR2, SSR1, SS12, rho, n1*n2, *lambda_eps2, *lambda_eps1, a_eps, b_eps);
  
}
