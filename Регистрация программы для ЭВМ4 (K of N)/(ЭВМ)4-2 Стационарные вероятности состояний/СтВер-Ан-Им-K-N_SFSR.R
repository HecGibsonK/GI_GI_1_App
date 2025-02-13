########### <GI_K/GI/1> ##################


rm(list=ls())
require(graphics)
require(xtable)
require(tcltk)
require(VGAM)


###############################################################################################
###SIMULATION MODEL (ИМИТАЦИОННАЯ МОДЕЛЬ) #####################################################
###############################################################################################

Sys_process=function(a1, b1, K, N, TT, NG,SystemType=NULL){
  
  #a1	- Average time between element failures (Среднее время между отказами элементов)
  #b1	- Average repair time (Среднее время ремонта)
  #K	- The number of elements operating in the system (Число элементов обеспечивающих оперативное функционирование системы)
  #N	- The number of elements in the system (Число элементов в системе)
  #TT	- Maximum model run time (Максимальное модельное время прогона)
  #NG	- Number of Trajectory Graphics (Число Графиков траекторий)
  
  alfa=1/a1; #The failure rate between element failures (Интенсивность отказа между отказами элементов)
  
  beta=1/b1; #Intensity of repair of the exponential distribution (Интенсивносность ремонта экспоненциального распределения)
  
  
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
    
    #cat('времени отказа элементов s: \n')
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
      
      #cat('Генерация таблицы времени в системе E0: \n')
      #print(c(E0))
      #cat('Генерация вектор времени в системе tv: \n')
      #print(tv)
      #cat('Генерация времени затрачно в системе tp: \n')
      #print(tp)
      #cat('времени отказа элементов в порядке snv : \n')
      #print(snv)
      #cat('Генерация времени ремонта sr: \n')
      #print(sr)
      
      r=array(r, dim=c(1,3,k));
      r[,,k] <-c(t,i,j); 
      i=j;
      k=k+1;
    }
  }
  
  ##
  if(SystemType=='<Expo_k/Weibull(W)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rexp(N,alfa);
    
    #cat('времени отказа элементов s: \n')
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
            if(i==1){sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rexp(1,alfa);
            snv = c(s_r[s_r>0], s_new);
            sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));
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
      
      #cat('Генерация таблицы времени в системе E0: \n')
      #print(c(E0))
      #cat('Генерация вектор времени в системе tv: \n')
      #print(tv)
      #cat('Генерация времени затрачно в системе tp: \n')
      #print(tp)
      #cat('времени отказа элементов в порядке snv : \n')
      #print(snv)
      #cat('Генерация времени ремонта sr: \n')
      #print(sr)
      
      r=array(r, dim=c(1,3,k));
      r[,,k] <-c(t,i,j); 
      i=j;
      k=k+1;
    }
  }
  
  if(SystemType=='<Expo_k/Pareto(P)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rexp(N,alfa);
    
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
            if(i==1){sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rexp(1,alfa);
            snv = c(s_r[s_r>0], s_new);
            sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);
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
  
  if(SystemType=='<Expo_k/Gamma(G)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rexp(N,alfa);
    
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
          sr = 0;
          s_r = sort(snv[] - tp);
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
            if(i==1){sr=rgamma(1,shape=G,scale = b1/G);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rexp(1,alfa);
            snv = c(s_r[s_r>0], s_new);
            sr=rgamma(1,shape=G,scale = b1/G);
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
  
  if(SystemType=='<Expo_k/Lognormal(sig)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rexp(N,alfa);
    b=log(b1)-(sig^2)/2;
    
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
            if(i==1){sr=rlnorm(1, meanlog =b, sdlog = sig);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rexp(1,alfa);
            snv = c(s_r[s_r>0], s_new);
            sr=rlnorm(1, meanlog =b, sdlog = sig);
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
  if(SystemType=='<Weibull(W)_k/Expo/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
    
    while(t<TT){
      
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
          s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
            s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
  
  if(SystemType=='<Weibull(W)_k/Weibull(W)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
    
    while(t<TT){
      
      if(i==0){
        tv = append(tv, t);
        if(length(tv)==1){
          tp = 0; sr = 0; snv=sort(s);
          sn=min(snv); t_otk=t+sn;
          t_rem=Inf; j=i+1; t=t_otk;
        }
        else{
          tp = tv[length(tv)] - tv[length(tv)-1];
          sr = 0; s_r = sort(snv[] - tp);
          s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
            if(i==1){sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
            snv = c(s_r[s_r>0], s_new);
            sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));
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
  
  if(SystemType=='<Weibull(W)_k/Pareto(P)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
    
    while(t<TT){
      
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
          s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
            if(i==1){sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
            snv = c(s_r[s_r>0], s_new);
            sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);
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
  
  if(SystemType=='<Weibull(W)_k/Gamma(G)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
    
    while(t<TT){
      
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
          s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
            if(i==1){sr=rgamma(1,shape=G,scale = b1/G);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
            snv = c(s_r[s_r>0], s_new);
            sr=rgamma(1,shape=G,scale = b1/G);
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
  
  if(SystemType=='<Weibull(W)_k/Lognormal(sig)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
    b=log(b1)-(sig^2)/2;
    
    while(t<TT){
      
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
          s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
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
            if(i==1){sr=rlnorm(1, meanlog =b, sdlog = sig);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
            snv = c(s_r[s_r>0], s_new);
            sr=rlnorm(1, meanlog =b, sdlog = sig);
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
  if(SystemType=='<Pareto(P)_k/Expo/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rpareto(N,shape=P,scale=(a1*(P-1))/P);
    
    while(t<TT){
      
      if(i==0){
        tv = append(tv, t);
        if(length(tv)==1){
          tp = 0; sr = 0; snv=sort(s);
          sn=min(snv);t_otk=t+sn;
          t_rem=Inf; j=i+1; t=t_otk;
        }
        else{
          tp = tv[length(tv)] - tv[length(tv)-1];
          sr = 0;	s_r = sort(snv[] - tp);
          s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
            s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
  
  if(SystemType=='<Pareto(P)_k/Weibull(W)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rpareto(N,shape=P,scale=(a1*(P-1))/P);
    
    while(t<TT){
      
      if(i==0){
        tv = append(tv, t);
        if(length(tv)==1){
          tp = 0; sr = 0; snv=sort(s);
          sn=min(snv);t_otk=t+sn;
          t_rem=Inf; j=i+1; t=t_otk;
        }
        else{
          tp = tv[length(tv)] - tv[length(tv)-1];
          sr = 0;	s_r = sort(snv[] - tp);
          s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
            if(i==1){sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
            snv = c(s_r[s_r>0], s_new);
            sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));
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
  
  if(SystemType=='<Pareto(P)_k/Pareto(P)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rpareto(N,shape=P,scale=(a1*(P-1))/P);
    
    while(t<TT){
      
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
          s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
            if(i==1){sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
            snv = c(s_r[s_r>0], s_new);
            sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);
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
  
  if(SystemType=='<Pareto(P)_k/Gamma(G)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rpareto(N,shape=P,scale=(a1*(P-1))/P);
    
    while(t<TT){
      
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
          s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
            if(i==1){sr=rgamma(1,shape=G,scale = b1/G);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
            snv = c(s_r[s_r>0], s_new);
            sr=rgamma(1,shape=G,scale = b1/G);
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
  
  if(SystemType=='<Pareto(P)_k/Lognormal(sig)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rpareto(N,shape=P,scale=(a1*(P-1))/P);
    b=log(b1)-(sig^2)/2;
    
    while(t<TT){
      
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
          s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
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
            if(i==1){sr=rlnorm(1, meanlog =b, sdlog = sig);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
            snv = c(s_r[s_r>0], s_new);
            sr=rlnorm(1, meanlog =b, sdlog = sig);
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
  if(SystemType=='<Gamma(G)_k/Expo/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rgamma(N,shape=G,scale = a1/G);
    
    while(t<TT){
      
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
          s_new = rgamma(1,shape=G,scale = a1/G);
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
            s_new = rgamma(1,shape=G,scale = a1/G);
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
  
  if(SystemType=='<Gamma(G)_k/Weibull(W)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rgamma(N,shape=G,scale = a1/G);
    
    while(t<TT){
      
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
          s_new = rgamma(1,shape=G,scale = a1/G);
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
            if(i==1){sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rgamma(1,shape=G,scale = a1/G);
            snv = c(s_r[s_r>0], s_new);
            sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));
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
  
  if(SystemType=='<Gamma(G)_k/Pareto(P)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rgamma(N,shape=G,scale = a1/G);
    
    while(t<TT){
      
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
          s_new = rgamma(1,shape=G,scale = a1/G);
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
            if(i==1){sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rgamma(1,shape=G,scale = a1/G);
            snv = c(s_r[s_r>0], s_new);
            sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);
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
  
  if(SystemType=='<Gamma(G)_k/Gamma(G)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rgamma(N,shape=G,scale = a1/G);
    
    while(t<TT){
      
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
          s_new = rgamma(1,shape=G,scale = a1/G);
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
            if(i==1){sr=rgamma(1,shape=G,scale = b1/G);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rgamma(1,shape=G,scale = a1/G);
            snv = c(s_r[s_r>0], s_new);
            sr=rgamma(1,shape=G,scale = b1/G);
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
  
  if(SystemType=='<Gamma(G)_k/Lognormal(sig)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    s=rgamma(N,shape=G,scale = a1/G);
    b=log(b1)-(sig^2)/2;
    
    while(t<TT){
      
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
          s_new = rgamma(1,shape=G,scale = a1/G);
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
            if(i==1){sr=rlnorm(1, meanlog =b, sdlog = sig);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rgamma(1,shape=G,scale = a1/G);
            snv = c(s_r[s_r>0], s_new);
            sr=rlnorm(1, meanlog =b, sdlog = sig);
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
  if(SystemType=='<Lognormal(sig)_k/Expo/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    a=log(a1)-(sig^2)/2;
    s=rlnorm(N, meanlog =a, sdlog = sig);
    
    while(t<TT){
      
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
          s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
            s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
  
  if(SystemType=='<Lognormal(sig)_k/Weibull(W)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    a=log(a1)-(sig^2)/2;
    s=rlnorm(N, meanlog =a, sdlog = sig);
    
    while(t<TT){
      
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
          s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
            if(i==1){sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rlnorm(1, meanlog =a, sdlog = sig);
            snv = c(s_r[s_r>0], s_new);
            sr=rweibull(1,shape=W,scale=b1/gamma(1+(1/W)));
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
  
  if(SystemType=='<Lognormal(sig)_k/Pareto(P)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    a=log(a1)-(sig^2)/2;
    s=rlnorm(N, meanlog =a, sdlog = sig);
    
    while(t<TT){
      
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
          s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
            if(i==1){sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rlnorm(1, meanlog =a, sdlog = sig);
            snv = c(s_r[s_r>0], s_new);
            sr=rpareto(1,shape=P,scale=(b1*(P-1))/P);
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
  
  if(SystemType=='<Lognormal(sig)_k/Gamma(G)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    a=log(a1)-(sig^2)/2;
    s=rlnorm(N, meanlog =a, sdlog = sig);
    
    while(t<TT){
      
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
          s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
            if(i==1){sr=rgamma(1,shape=G,scale = b1/G);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rlnorm(1, meanlog =a, sdlog = sig);
            snv = c(s_r[s_r>0], s_new);
            sr=rgamma(1,shape=G,scale = b1/G);
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
  
  if(SystemType=='<Lognormal(sig)_k/Lognormal(sig)/1>'){
    s=vector('numeric',N);
    snv=vector('numeric',N);
    tv=vector();
    sr=vector();
    
    a=log(a1)-(sig^2)/2;
    s=rlnorm(N, meanlog =a, sdlog = sig);
    b=log(b1)-(sig_b^2)/2;
    
    while(t<TT){
      
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
          s_new = rlnorm(1, meanlog =a, sdlog = sig);
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
            if(i==1){sr=rlnorm(1, meanlog =b, sdlog = sig);}
            else if(length(tv)>2){
              sr = abs(sr[] - tp); 
            }
          }
          
          if(t==t_rem){
            s_new = rlnorm(1, meanlog =a, sdlog = sig);
            snv = c(s_r[s_r>0], s_new);
            sr=rlnorm(1, meanlog =b, sdlog = sig);
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
####ANALYTICAL MODEL (АНАЛИТИЧЕСКАЯ МОДЕЛЬ) ####################################################
################################################################################################

SysA_process=function(SystemType=NULL){
  
  L=vector('numeric',K);
  
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){
    a=1/a1;
    for(i in 0:(K-1)){
      L[i+1] = a*(N-i); # интенсивновность отказа элементов Экспоненциального распределения облег. рез.
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
  
  if(SystemType=='<Expo_k/Weibull(W)/1>'){
    P_WB = vector('numeric',K);
    A = vector('numeric',K);
    C1 = vector('numeric',K+1);
    Pi_WB = vector('numeric',K+1);
    
    b=b1/gamma(1+(1/W));
    
    F1 = function(x){dweibull(x, shape = W, scale = b, log = FALSE)} 
    #F1 = function(x){(W/b)*((x/b)^(W-1))*exp(-(x/b)^W)}
    
    for(i in 0:(K-1)){
      P_WB[i+1] = integrate(function(x){exp(-L[i+1]*x)*F1(x)},0,Inf)$value
    }
    
    #b1 = integrate(function(x){x*F1(x)},0,Inf)$value
    #b1 = b*gamma(1+(1/W));
    
    if(K==2){
      A[1] = 1;
      A[K]= A[K-1];
    }
    
    if(K>=3){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_WB[2])*(1/P_WB[3]);
      
      AK1=0; AK2=0;
      for(j in 1:(K-2)){
        AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
        AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_WB[j+1]);	
      }
      A[K] = A[K-1]*(1+P_WB[K]) + AK1 - AK2;
    }
    
    if(K>=4){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_WB[2])*(1/P_WB[3]);
      
      for(i in 2:(K-2)){
        AI1=0;
        for(j in 1:(i-1)){
          AI1 =AI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]);
        }
        AI2=0;
        for(j in 1:i){
          AI2 =AI2 + sum(((-1)^(i+1-j))*(prod(L[(j:i)+1]/(L[j+1]-L[(j:i)+2])))*A[j]*P_WB[j+1]);
        }
        A[i+1] = (A[i] + AI1 - AI2)*(1/P_WB[i+2]);
        AK1=0; AK2=0;
        for(j in 1:(K-2)){
          AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
          AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_WB[j+1]);	
        }
      }
      A[K] = A[K-1]*(1+P_WB[K]) + AK1 - AK2;
    }
    
    C1[1] = P_WB[2]/L[1];
    C1[2] = (1 - P_WB[2])/L[2];
    
    if(K>=3){
      for(i in 2:(K-1)){
        CI1=0;
        for(j in 1:(i-1)){
          CI1 =CI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_WB[j+1])/L[j+1]));	
        }
        C1[i+1] = A[i]*((1-P_WB[i+1])/L[i+1]) + CI1;
      }
    }
    Cn=0;
    if(K==2){
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_WB[K])/L[K]);
    }
    else{
      for(j in 1:(K-2)){
        Cn =Cn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_WB[j+1])/L[j+1]));
      }
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_WB[K])/L[K]) + Cn;
    }
    
    C = (sum(C1))^-1;
    
    Pi_WB[1] = C*(P_WB[2]/L[1]);
    Pi_WB[2] = C*((1 - P_WB[2])/L[2]);
    
    if(K>2){
      for(i in 2:(K-1)){
        Pi=0;
        for(j in 1:(i-1)){
          Pi =Pi + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_WB[j+1])/L[j+1]));
        }
        Pi_WB[i+1] = C*(A[i]*((1-P_WB[i+1])/L[i+1]) + Pi);
      }
    }
    
    Pn=0;
    if(K==2){
      Pi_WB[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_WB[K])/L[K]));
    }
    else{ 
      for(j in 1:(K-2)){
        Pn =Pn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_WB[j+1])/L[j+1]));
      }
      Pi_WB[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_WB[K])/L[K]) + Pn);
    }
    Pi = Pi_WB;
  }
  
  if(SystemType=='<Expo_k/Pareto(P)/1>'){
    P_PAR = vector('numeric',K);
    A = vector('numeric',K);
    C1 = vector('numeric',K+1);
    Pi_PAR = vector('numeric',K+1);
    
    b=(b1*(P-1))/P; #### P>1; else "Inf"
    F1 = function(x){dpareto(x, shape = P, scale = b, log = FALSE)} 
    #F1 = function(x){(P*(b^P))/(x^(P+1))}
    
    for(i in 0:(K-1)){
      P_PAR[i+1] = integrate(function(x){exp(-L[i+1]*x)*F1(x)},b,Inf)$value
    }
    
    #b1 = integrate(function(x){x*F1(x)},0,Inf)$value
    #b1 = (P*b)/(P-1); #### P>1; else "Inf"
    
    if(K==2){
      A[1] = 1;
      A[K]= A[K-1];
    }
    
    if(K>=3){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_PAR[2])*(1/P_PAR[3]);
      
      AK1=0; AK2=0;
      for(j in 1:(K-2)){
        AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
        AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_PAR[j+1]);	
      }
      A[K] = A[K-1]*(1+P_PAR[K]) + AK1 - AK2;
    }
    
    if(K>=4){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_PAR[2])*(1/P_PAR[3]);
      
      for(i in 2:(K-2)){
        AI1=0;
        for(j in 1:(i-1)){
          AI1 =AI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]);
        }
        AI2=0;
        for(j in 1:i){
          AI2 =AI2 + sum(((-1)^(i+1-j))*(prod(L[(j:i)+1]/(L[j+1]-L[(j:i)+2])))*A[j]*P_PAR[j+1]);
        }
        A[i+1] = (A[i] + AI1 - AI2)*(1/P_PAR[i+2]);
        AK1=0; AK2=0;
        for(j in 1:(K-2)){
          AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
          AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_PAR[j+1]);	
        }
      }
      A[K] = A[K-1]*(1+P_PAR[K]) + AK1 - AK2;
    }
    
    C1[1] = P_PAR[2]/L[1];
    C1[2] = (1 - P_PAR[2])/L[2];
    
    if(K>=3){
      for(i in 2:(K-1)){
        CI1=0;
        for(j in 1:(i-1)){
          CI1 =CI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_PAR[j+1])/L[j+1]));	
        }
        C1[i+1] = A[i]*((1-P_PAR[i+1])/L[i+1]) + CI1;
      }
    }
    Cn=0;
    if(K==2){
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_PAR[K])/L[K]);
    }
    else{
      for(j in 1:(K-2)){
        Cn =Cn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_PAR[j+1])/L[j+1]));
      }
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_PAR[K])/L[K]) + Cn;
    }
    
    C = (sum(C1))^-1;
    
    Pi_PAR[1] = C*(P_PAR[2]/L[1]);
    Pi_PAR[2] = C*((1 - P_PAR[2])/L[2]);
    
    if(K>2){
      for(i in 2:(K-1)){
        Pi=0;
        for(j in 1:(i-1)){
          Pi =Pi + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_PAR[j+1])/L[j+1]));
        }
        Pi_PAR[i+1] = C*(A[i]*((1-P_PAR[i+1])/L[i+1]) + Pi);
      }
    }
    
    Pn=0;
    if(K==2){
      Pi_PAR[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_PAR[K])/L[K]));
    }
    else{ 
      for(j in 1:(K-2)){
        Pn =Pn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_PAR[j+1])/L[j+1]));
      }
      Pi_PAR[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_PAR[K])/L[K]) + Pn);
    }
    
    Pi = Pi_PAR;
  }
  
  if(SystemType=='<Expo_k/Gamma(G)/1>'){
    P_G = vector('numeric',K);
    A = vector('numeric',K);
    C1 = vector('numeric',K+1);
    Pi_G = vector('numeric',K+1);
    
    b=b1/G;
    
    F1 = function(x){dgamma(x, shape = G, scale = b, log = FALSE)} 
    #F1 = function(x){(1/(gamma(G)*b^G))*(x^(G-1) * exp(-x/b))}
    
    for(i in 0:(K-1)){
      P_G[i+1] = integrate(function(x){exp(-L[i+1]*x)*F1(x)},0,Inf)$value
    }
    
    #b1 = integrate(function(x){x*F1(x)},0,Inf)$value
    #b1 = b*G
    
    if(K==2){
      A[1] = 1;
      A[K]= A[K-1];
    }
    
    if(K>=3){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_G[2])*(1/P_G[3]);
      
      AK1=0; AK2=0;
      for(j in 1:(K-2)){
        AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
        AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_G[j+1]);	
      }
      A[K] = A[K-1]*(1+P_G[K]) + AK1 - AK2;
    }
    
    if(K>=4){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_G[2])*(1/P_G[3]);
      
      for(i in 2:(K-2)){
        AI1=0;
        for(j in 1:(i-1)){
          AI1 =AI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]);
        }
        AI2=0;
        for(j in 1:i){
          AI2 =AI2 + sum(((-1)^(i+1-j))*(prod(L[(j:i)+1]/(L[j+1]-L[(j:i)+2])))*A[j]*P_G[j+1]);
        }
        A[i+1] = (A[i] + AI1 - AI2)*(1/P_G[i+2]);
        AK1=0; AK2=0;
        for(j in 1:(K-2)){
          AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
          AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_G[j+1]);	
        }
      }
      A[K] = A[K-1]*(1+P_G[K]) + AK1 - AK2;
    }
    
    C1[1] = P_G[2]/L[1];
    C1[2] = (1 - P_G[2])/L[2];
    
    if(K>=3){
      for(i in 2:(K-1)){
        CI1=0;
        for(j in 1:(i-1)){
          CI1 =CI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_G[j+1])/L[j+1]));	
        }
        C1[i+1] = A[i]*((1-P_G[i+1])/L[i+1]) + CI1;
      }
    }
    Cn=0;
    if(K==2){
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_G[K])/L[K]);
    }
    else{
      for(j in 1:(K-2)){
        Cn =Cn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_G[j+1])/L[j+1]));
      }
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_G[K])/L[K]) + Cn;
    }
    
    C = (sum(C1))^-1;
    
    Pi_G[1] = C*(P_G[2]/L[1]);
    Pi_G[2] = C*((1 - P_G[2])/L[2]);
    
    if(K>2){
      for(i in 2:(K-1)){
        Pi=0;
        for(j in 1:(i-1)){
          Pi =Pi + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_G[j+1])/L[j+1]));
        }
        Pi_G[i+1] = C*(A[i]*((1-P_G[i+1])/L[i+1]) + Pi);
      }
    }
    
    Pn=0;
    if(K==2){
      Pi_G[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_G[K])/L[K]));
    }
    else{ 
      for(j in 1:(K-2)){
        Pn =Pn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_G[j+1])/L[j+1]));
      }
      Pi_G[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_G[K])/L[K]) + Pn);
    }
    
    Pi = Pi_G;	
  }
  
  if(SystemType=='<Expo_k/Lognormal(sig)/1>'){
    P_LN = vector('numeric',K);
    A = vector('numeric',K);
    C1 = vector('numeric',K+1);
    Pi_LN = vector('numeric',K+1);
    
    b=log(b1)-(sig^2)/2;
    
    F1 = function(x){dlnorm(x, meanlog =b, sdlog = sig, log = FALSE)} 
    #F1 = function(x){(1/(x*sig*sqrt(2*pi)))*exp(-((log(x)-b)^2 / (2*sig^2)))}
    
    for(i in 0:(K-1)){
      P_LN[i+1] = integrate(function(x){exp(-L[i+1]*x)*F1(x)},0,Inf)$value
    }
    
    #b1 = integrate(function(x){x*F1(x)},0,Inf)$value
    #b1 = exp(b+(sig^2/2))
    
    if(K==2){
      A[1] = 1;
      A[K]= A[K-1];
    }
    
    if(K>=3){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_LN[2])*(1/P_LN[3]);
      
      AK1=0; AK2=0;
      for(j in 1:(K-2)){
        AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
        AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_LN[j+1]);	
      }
      A[K] = A[K-1]*(1+P_LN[K]) + AK1 - AK2;
    }
    
    if(K>=4){
      A[1] = 1;
      A[2] = (1 - (1-(L[2]/(L[2]-L[3])))*P_LN[2])*(1/P_LN[3]);
      
      for(i in 2:(K-2)){
        AI1=0;
        for(j in 1:(i-1)){
          AI1 =AI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]);
        }
        AI2=0;
        for(j in 1:i){
          AI2 =AI2 + sum(((-1)^(i+1-j))*(prod(L[(j:i)+1]/(L[j+1]-L[(j:i)+2])))*A[j]*P_LN[j+1]);
        }
        A[i+1] = (A[i] + AI1 - AI2)*(1/P_LN[i+2]);
        AK1=0; AK2=0;
        for(j in 1:(K-2)){
          AK1 =AK1 + sum(((-1)^(K-1-j))*(prod(L[(j:(K-2))+1]/(L[j+1]-L[(j:(K-2))+2])))*A[j]);
          AK2 =AK2 + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*P_LN[j+1]);	
        }
      }
      A[K] = A[K-1]*(1+P_LN[K]) + AK1 - AK2;
    }
    
    C1[1] = P_LN[2]/L[1];
    C1[2] = (1 - P_LN[2])/L[2];
    
    if(K>=3){
      for(i in 2:(K-1)){
        CI1=0;
        for(j in 1:(i-1)){
          CI1 =CI1 + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_LN[j+1])/L[j+1]));	
        }
        C1[i+1] = A[i]*((1-P_LN[i+1])/L[i+1]) + CI1;
      }
    }
    Cn=0;
    if(K==2){
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_LN[K])/L[K]);
    }
    else{
      for(j in 1:(K-2)){
        Cn =Cn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_LN[j+1])/L[j+1]));
      }
      C1[K+1] = A[K]*b1 - A[K-1]*((1-P_LN[K])/L[K]) + Cn;
    }
    
    C = (sum(C1))^-1;
    
    Pi_LN[1] = C*(P_LN[2]/L[1]);
    Pi_LN[2] = C*((1 - P_LN[2])/L[2]);
    
    if(K>2){
      for(i in 2:(K-1)){
        Pi=0;
        for(j in 1:(i-1)){
          Pi =Pi + sum(((-1)^(i-j))*(prod(L[(j:(i-1))+1]/(L[j+1]-L[(j:(i-1))+2])))*A[j]*((1-P_LN[j+1])/L[j+1]));
        }
        Pi_LN[i+1] = C*(A[i]*((1-P_LN[i+1])/L[i+1]) + Pi);
      }
    }
    
    Pn=0;
    if(K==2){
      Pi_LN[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_LN[K])/L[K]));
    }
    else{ 
      for(j in 1:(K-2)){
        Pn =Pn + sum(((-1)^(K-j))*(prod(L[(j:(K-2))+2]/(L[j+1]-L[(j:(K-2))+2])))*A[j]*((1-P_LN[j+1])/L[j+1]));
      }
      Pi_LN[K+1] = C*(A[K]*b1 - A[K-1]*((1-P_LN[K])/L[K]) + Pn);
    }
    Pi = Pi_LN;	
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
    cat('input parameter (вводите параметр) a1 "Average time between element failures (Среднее время между отказами элементов)": \n');a1<<-scan(nlines=1);
    cat('input parameter (вводите параметр) b1 "Average repair time (Среднее время ремонта)": \n');b1<<-scan(nlines=1);
    cat('input parameter (вводите параметр) K "The number of elements operating in the system (Число элементов действующие в системе)": \n');K<<-scan(nlines=1);
    cat('input parameter (вводите параметр) N "The number of elements in the system (Число элементов в системе)": \n');N<<-scan(nlines=1);
    cat('input parameter (вводите параметр) TT "Maximum model run time (Максимальное модельное время прогона)": \n');TT<<-scan(nlines=1);
    cat('input parameter (вводите параметр) NG "Number of Trajectory Graphics (Число Графиков траекторий)": \n');NG<<-scan(nlines=1);
  }
  
  if(SystemType=='<Expo_k/Weibull(W)/1>'){
    cat('input parameter (вводите параметр) W "Parameter repair time distribution element (параметр распределения времени ремонта элемента)": \n');W<<-scan(nlines=1);
  }
  if(SystemType=='<Weibull(W)_k/Expo/1>'){
    cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1);
  }
  if(SystemType=='<Expo_k/Pareto(P)/1>'){
    cat('input parameter (вводите параметр) P "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');P<<-scan(nlines=1);
  }
  if(SystemType=='<Pareto(P)_k/Expo/1>'){
    cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1);
  }
  if(SystemType=='<Expo_k/Gamma(G)/1>'){
    cat('input parameter (вводите параметр) G "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');G<<-scan(nlines=1)
  }
  if(SystemType=='<Gamma(G)_k/Expo/1>'){
    cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1)
  }
  if(SystemType=='<Expo_k/Lognormal(sig)/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');sig<<-scan(nlines=1);
  }
  if(SystemType=='<Lognormal(sig)_k/Expo/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1);
  }
  #	
  if(SystemType=='<Weibull(W)_k/Weibull(W)/1>'){
    cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) W_b "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');W_b<<-scan(nlines=1);
  }
  if(SystemType=='<Weibull(W)_k/Pareto(P)/1>'){
    cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) P "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');P<<-scan(nlines=1);
  }
  if(SystemType=='<Weibull(W)_k/Gamma(G)/1>'){
    cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) G "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');G<<-scan(nlines=1);
  }
  if(SystemType=='<Weibull(W)_k/Lognormal(sig)/1>'){
    cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) sig "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');sig<<-scan(nlines=1);
  }
  #
  if(SystemType=='<Pareto(P)_k/Weibull(W)/1>'){
    
    cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) W "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');W<<-scan(nlines=1);
  }
  if(SystemType=='<Pareto(P)_k/Pareto(P)/1>'){
    cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1);
    cat('input parameter (вводите параметр) P_b "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');P_b<<-scan(nlines=1);
  }
  if(SystemType=='<Pareto(P)_k/Gamma(G)/1>'){
    cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1);
    cat('input parameter (вводите параметр) G "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');G<<-scan(nlines=1);
  }
  if(SystemType=='<Pareto(P)_k/Lognormal(sig)/1>'){
    cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) sig "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');sig<<-scan(nlines=1);
  }
  #
  if(SystemType=='<Gamma(G)_k/Weibull(W)/1>'){
    cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) W "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');W<<-scan(nlines=1);
  }
  if(SystemType=='<Gamma(G)_k/Pareto(P)/1>'){
    cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) P "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');P<<-scan(nlines=1);
  }
  if(SystemType=='<Gamma(G)_k/Gamma(G)/1>'){
    cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) G_b "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');G_b<<-scan(nlines=1);
  }
  if(SystemType=='<Gamma(G)_k/Lognormal(sig)/1>'){
    cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) sig "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');sig<<-scan(nlines=1);
  }
  #
  if(SystemType=='<Lognormal(sig)_k/Weibull(W)/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) W "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');W<<-scan(nlines=1);
  }
  if(SystemType=='<Lognormal(sig)_k/Pareto(P)/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) P "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');P<<-scan(nlines=1);
  }
  if(SystemType=='<Lognormal(sig)_k/Gamma(G)/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) G "repair time distribution parameter (параметр распределения времени ремонта)": \n');G<<-scan(nlines=1);
  }
  if(SystemType=='<Lognormal(sig)_k/Lognormal(sig)/1>'){
    cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1); 
    cat('input parameter (вводите параметр) sig_b "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');sig_b<<-scan(nlines=1);
  }
  #
  if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){
    Pi_a = SysA_process(SystemType);
  }
  
  ###########################################################################
  ## Графики траекторий
  
  pr=Sys_process(a1, b1, K, N, TT, NG,SystemType);
  
  NumF=pr[1,3,]
  plot(c(0,pr[1,1,]),c(0,NumF),col="green",type="s", yaxt='n', xaxt='n',xlab="", ylab="X(t)")
  mtext("t",side=1,at=T+5,line=1,cex=1.1)
  axis(side = 2, at = seq(0,K), las=1)
  axis(side = 1, at = c(0,pr[1,1,]), las=2)
  abline(v=TT,col="red")
  
  ###########################################################################
  ## Стационарные Вероятности состояния системы
  
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
    Pi = matrix(c(PIV,Pi_a), nrow = 2, ncol = (K+1), byrow=T,list(c('Simul.(Имит.)','Analit.(Аналит.)'))) 
    
    cat('Allowable error (Допустимая погрешность): \n')
    PP = abs(PIV-Pi_a)
    #ss = max(PP);ss  
    print(matrix(c(PP), nrow = 1, ncol = (K+1)))
  }
  else{
    Pi = matrix(c(PIV), nrow = 1, ncol = (K+1), byrow=T,list(c('Simul.(Имит.)')))
  }
  
  return(Pi)
}

##############################################################################################
##############################################################################################

#Pi = main();Pi

