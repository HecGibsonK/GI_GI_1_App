########### <GI_K/GI/1> ##################


rm(list=ls())
require(graphics)
require(xtable)
require(tcltk)
require(VGAM)


###############################################################################################
###SIMULATION MODEL (???????????? ??????) #####################################################
###############################################################################################

Sys_process=function(a1, b1, K, N, TT, NG,SystemType=NULL){
  
  #a1	- Average time between element failures (??????? ????? ????? ???????? ?????????)
  #b1	- Average repair time (??????? ????? ???????)
  #K	- The number of elements operating in the system (????? ????????? ?????????????? ??????????? ???????????????? ???????)
  #N	- The number of elements in the system (????? ????????? ? ???????)
  #TT	- Maximum model run time (???????????? ????????? ????? ???????)
  #NG	- Number of Trajectory Graphics (????? ???????? ??????????)
  
  alfa=1/a1; #The failure rate between element failures (????????????? ?????? ????? ???????? ?????????)
  
  beta=1/b1; #Intensity of repair of the exponential distribution (???????????????? ??????? ????????????????? ?????????????)
  
  
  r=array(, dim=c(1,3,N));
  
  i=0;
  t=0;
  tv=0; tv = append(tv, tv);
  t_otk=0; t_rem=Inf;
  k=1;
  
  ##
  if(SystemType=='<Expo_k/Expo/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rexp(N,alfa);
    
    #cat('??????? ?????? ????????? s: \n')
    #print(s)
    
    while(t<TT){
      
      E0=array(, dim=c(1,1,1));
      nl=1;
      E0=array(E0, dim=c(1,1,nl));
      E0[,,nl]=c(t);
      
      nl=nl+1;
      
      if(i==0){
        tv = append(tv, t);
        if(length(tv)==1){
          tp = 0; sr = 0; snv=sort(s);
          sn=min(snv);t_otk=t+sn;
          t_rem=Inf; j=i+1; t=t_otk;
        }
        else{
          tp = tv[length(tv)] - tv[length(tv)-1];
          sr = 0; s_r = sort(snv[] - tp);
          s_new = rexp(1,alfa);
          snv = c(s_r[s_r>0], s_new);
          sn=min(snv); t_otk=t+sn;
          t_rem=Inf; j=i+1; t=t_otk;
        }
      }
      
      for(v in 1:(K-1)){ 
        if(i==v){
          tv = append(tv, t);					
          tp = tv[length(tv)] - tv[length(tv)-1];
          s_r = sort(snv[] - tp);
          
          if(t==t_otk){
            snv = s_r[-1];
            if(i==1){sr=rexp(1,beta);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rexp(1,alfa);
            snv = c(s_r[s_r>0], s_new);
            sr=rexp(1,beta);
          }
          
          sn=min(snv); 
          t_otk=t+sn; t_rem=t+sr; 
          if(t_otk<t_rem){
            j=i+1; t=t_otk;
          }
          else{
            j=i-1; t=t_rem;
          }
        }
      }
      if(i==K){
        tv = append(tv, t);
        tp = tv[length(tv)] - tv[length(tv)-1];
        
        if(K==N){
          snv = 0;
          sr = abs(sr[] - tp); t_rem=t+sr;
        }
        else{
          s_r = sort(snv[] - tp); sr = abs(sr[] - tp); 
          t_rem = t+sr; snv = s_r[-1] + sr;
        }
        t_otk = Inf; j = i-1; t = t_rem;
      }
      
      if(t>TT){t=TT}
      
      r=array(r, dim=c(1,3,k));
      r[,,k] <-c(t,i,j); 
      i=j;
      k=k+1;
    }
  }
  
  ##
  return(r)
}

################################################################################################
####ANALYTICAL MODEL (????????????? ??????) ####################################################
################################################################################################

SysA_process=function(SystemType=NULL){
  
  L=vector('numeric',K);
  
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){
    a=1/a1;
    for(i in 0:(K-1)){
      L[i+1] = a*(N-i); # ???????????????? ?????? ????????? ????????????????? ????????????? ?????. ???.
    }
  }
  
  if(SystemType=='<Expo_k/Expo/1>'){
    P_EXP = vector('numeric',K);
    A = vector('numeric',K);
    C1 = vector('numeric',K+1);
    Pi_EXP = vector('numeric',K+1);
    
    b=1/b1;
    
    F1 = function(x){dexp(x, rate = b, log = FALSE)} 
    #F1 = function(x){b*exp(-b*x)} 
    
    for(i in 0:(K-1)){
      P_EXP[i+1] = integrate(function(x){exp(-L[i+1]*x)*F1(x)},0,Inf)$value
    }
    
    #b1 = integrate(function(x){x*F1(x)},0,Inf)$value
    #b1 = 1/b
    
    if(K==2){
      A[1] = 1;
      A[K]= A[K-1];
    }
    
    if(K>=3){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_EXP[2])*(1/P_EXP[3]);
      
      AK1=0; AK2=0;
      for(j in 1:(K-2)){
        AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
        AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_EXP[j+1]);	
      }
      A[K] = A[K-1]*(1+P_EXP[K]) + AK1 - AK2;
    }
    
    if(K>=4){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_EXP[2])*(1/P_EXP[3]);
      
      for(i in 2:(K-2)){
        AI1=0;
        for(j in 1:(i-1)){
          AI1 =AI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]);
        }
        AI2=0;
        for(j in 1:i){
          AI2 =AI2 + sum(((-1)^(i+1-j))*(prod(L[(j:i)+1]/(L[j+1]-L[(j:i)+2])))*A[j]*P_EXP[j+1]);
        }
        A[i+1] = (A[i] + AI1 - AI2)*(1/P_EXP[i+2]);
        AK1=0; AK2=0;
        for(j in 1:(K-2)){
          AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
          AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_EXP[j+1]);	
        }
      }
      A[K] = A[K-1]*(1+P_EXP[K]) + AK1 - AK2;
    }
    
    C1[1] = P_EXP[2]/L[1];
    C1[2] = (1 - P_EXP[2])/L[2];
    
    if(K>=3){
      for(i in 2:(K-1)){
        CI1=0;
        for(j in 1:(i-1)){
          CI1 =CI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_EXP[j+1])/L[j+1]));	
        }
        C1[i+1] = A[i]*((1-P_EXP[i+1])/L[i+1]) + CI1;
      }
    }
    Cn=0;
    if(K==2){
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_EXP[K])/L[K]);
    }
    else{
      for(j in 1:(K-2)){
        Cn =Cn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_EXP[j+1])/L[j+1]));
      }
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_EXP[K])/L[K]) + Cn;
    }
    
    C = (sum(C1))^-1;
    
    Pi_EXP[1] = C*(P_EXP[2]/L[1]);
    Pi_EXP[2] = C*((1 - P_EXP[2])/L[2]);
    
    if(K>2){
      for(i in 2:(K-1)){
        Pi=0;
        for(j in 1:(i-1)){
          Pi =Pi + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_EXP[j+1])/L[j+1]));
        }
        Pi_EXP[i+1] = C*(A[i]*((1-P_EXP[i+1])/L[i+1]) + Pi);
      }
    }
    
    Pn=0;
    if(K==2){
      Pi_EXP[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_EXP[K])/L[K]));
    }
    else{ 
      for(j in 1:(K-2)){
        Pn =Pn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_EXP[j+1])/L[j+1]));
      }
      Pi_EXP[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_EXP[K])/L[K]) + Pn);
    }
    Pi = Pi_EXP;	
  }
  
  return(Pi)
}

######################################################################################
######################################################################################

############install.packages("VGAM",dependencies=T,repos='http:// cran.rstudio.com/')
#### install.packages("installr"); require(installr); updateR();

main=function(SystemType=NULL){
  
  if(is.null(SystemType)==TRUE){
    SystemType=select.list(c('<Expo_k/Expo/1>','<Expo_k/Weibull(W)/1>','<Expo_k/Pareto(P)/1>','<Expo_k/Gamma(G)/1>','<Expo_k/Lognormal(sig)/1>',
                             '<Weibull(W)_k/Expo/1>','<Weibull(W)_k/Weibull(W)/1>','<Weibull(W)_k/Pareto(P)/1>','<Weibull(W)_k/Gamma(G)/1>','<Weibull(W)_k/Lognormal(sig)/1>',
                             '<Pareto(P)_k/Expo/1>','<Pareto(P)_k/Weibull(W)/1>','<Pareto(P)_k/Pareto(P)/1>','<Pareto(P)_k/Gamma(G)/1>','<Pareto(P)_k/Lognormal(sig)/1>',
                             '<Gamma(G)_k/Expo/1>','<Gamma(G)_k/Weibull(W)/1>','<Gamma(G)_k/Pareto(P)/1>','<Gamma(G)_k/Gamma(G)/1>','<Gamma(G)_k/Lognormal(sig)/1>',
                             '<Lognormal(sig)_k/Expo/1>','<Lognormal(sig)_k/Weibull(W)/1>','<Lognormal(sig)_k/Pareto(P)/1>','<Lognormal(sig)_k/Gamma(G)/1>','<Lognormal(sig)_k/Lognormal(sig)/1>'),
                           preselect = NULL, multiple = FALSE, title = NULL, graphics = getOption("menu.graphics"))
  }
  
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Weibull(W)_k/Expo/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Pareto(P)_k/Expo/1>'||
     SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Gamma(G)_k/Expo/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'||SystemType=='<Lognormal(sig)_k/Expo/1>'||
     SystemType=='<Weibull(W)_k/Weibull(W)/1>'||SystemType=='<Weibull(W)_k/Pareto(P)/1>'||SystemType=='<Weibull(W)_k/Gamma(G)/1>'||SystemType=='<Weibull(W)_k/Lognormal(sig)/1>'||
     SystemType=='<Pareto(P)_k/Weibull(W)/1>'||SystemType=='<Pareto(P)_k/Pareto(P)/1>'||SystemType=='<Pareto(P)_k/Gamma(G)/1>'||SystemType=='<Pareto(P)_k/Lognormal(sig)/1>'||
     SystemType=='<Gamma(G)_k/Weibull(W)/1>'||SystemType=='<Gamma(G)_k/Pareto(P)/1>'||SystemType=='<Gamma(G)_k/Gamma(G)/1>'||SystemType=='<Gamma(G)_k/Lognormal(sig)/1>'||
     SystemType=='<Lognormal(sig)_k/Weibull(W)/1>'||SystemType=='<Lognormal(sig)_k/Pareto(P)/1>'||SystemType=='<Lognormal(sig)_k/Gamma(G)/1>'||SystemType=='<Lognormal(sig)_k/Lognormal(sig)/1>'){
    cat('input parameter (??????? ????????) a1 "Average time between element failures (??????? ????? ????? ???????? ?????????)": \n');a1<<-scan(nlines=1);
    cat('input parameter (??????? ????????) b1 "Average repair time (??????? ????? ???????)": \n');b1<<-scan(nlines=1);
    cat('input parameter (??????? ????????) K "The number of elements operating in the system (????? ????????? ??????????? ? ???????)": \n');K<<-scan(nlines=1);
    cat('input parameter (??????? ????????) N "The number of elements in the system (????? ????????? ? ???????)": \n');N<<-scan(nlines=1);
    cat('input parameter (??????? ????????) TT "Maximum model run time (???????????? ????????? ????? ???????)": \n');TT<<-scan(nlines=1);
    cat('input parameter (??????? ????????) NG "Number of Trajectory Graphics (????? ???????? ??????????)": \n');NG<<-scan(nlines=1);
  }

  #
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){
    Pi_a = SysA_process(SystemType);
  }
  
  ###########################################################################
  ## ??????? ??????????
  
  pr=Sys_process(a1, b1, K, N, TT, NG,SystemType);
  
  NumF=pr[1,3,]
  plot(c(0,pr[1,1,]),c(0,NumF),col="green",type="s", yaxt='n', xaxt='n',xlab="", ylab="X(t)")
  mtext("t",side=1,at=T+5,line=1,cex=1.1)
  axis(side = 2, at = seq(0,K), las=1)
  axis(side = 1, at = c(0,pr[1,1,]), las=2)
  abline(v=TT,col="red")
  
  ###########################################################################
  ## ???????????? ??????????? ????????? ???????
  
  summ_T=function(pr,i){
    if(pr[1,2,1] == i){s= pr[1,1,1]}else{s= 0}
    if(dim(pr)[3]>1){
      for(k in 2:dim(pr)[3]){
        if(pr[1,2,k] == i){s= s+pr[1,1,k]-pr[1,1,k-1]}
      }
    }
    return(s)
  }
  
  viborka=vector();
  PIV=numeric()
  l=0;
  for(vv in 0:K){
    l=l+1
    for(j in 1:NG){
      
      prv=Sys_process(a1, b1, K, N, TT, NG,SystemType);
      
      pivv=vector('numeric',K+1)
      for(v in 0:(K)){
        pivv[v+1] = summ_T(prv,v)/TT;
      }
      
      viborka[j]=pivv[vv+1];
    }
    PIV[l] = mean(viborka); 
  }
  
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){
    Pi = matrix(c(PIV,Pi_a), nrow = 2, ncol = (K+1), byrow=T,list(c('Simul.(????.)','Analit.(??????.)'))) 
    
    cat('Allowable error (?????????? ???????????): \n')
    PP = abs(PIV-Pi_a)
    #ss = max(PP);ss  
    print(matrix(c(PP), nrow = 1, ncol = (K+1)))
  }
  else{
    Pi = matrix(c(PIV), nrow = 1, ncol = (K+1), byrow=T,list(c('Simul.(????.)')))
  }
  
  return(Pi)
}

##############################################################################################
##############################################################################################

#Pi = main();Pi

