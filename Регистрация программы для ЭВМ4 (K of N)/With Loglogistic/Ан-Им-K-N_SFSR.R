########### <GI_K/GI/1> ##################


rm(list=ls())
require(graphics)
require(xtable)
require(tcltk)
require(VGAM)
#require(fitdistrplus) # Just for Loglogistic regression using "dllogis"
#require(actuar) 		# Just for Loglogistic regression using "dllogis" 

#############################################
###SIMULATION MODEL (ИМИТАЦИОННАЯ МОДЕЛЬ) ###
#############################################

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
	
	if(SystemType=='<Expo_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		s=rexp(N,alfa);
		b = (b1*LL*sin(pi/LL))/pi;
		
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
						if(i==1){sr=rllogis(1, shape = LL, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rexp(1,alfa);
						snv = c(s_r[s_r>0], s_new);
						sr=rllogis(1, shape = LL, scale = b);
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
	
	if(SystemType=='<Weibull(W)_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		s=rweibull(N,shape=W,scale=a1/gamma(1+(1/W)));
		b = (b1*LL*sin(pi/LL))/pi;
		
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
						if(i==1){sr=rllogis(1, shape = LL, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rweibull(1,shape=W,scale=a1/gamma(1+(1/W)));
						snv = c(s_r[s_r>0], s_new);
						sr=rllogis(1, shape = LL, scale = b);
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
	
	if(SystemType=='<Pareto(P)_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		s = rpareto(N,shape=P,scale=(a1*(P-1))/P);
		b = (b1*LL*sin(pi/LL))/pi;
		
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
						if(i==1){sr=rllogis(1, shape = LL, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rpareto(1,shape=P,scale=(a1*(P-1))/P);
						snv = c(s_r[s_r>0], s_new);
						sr=rllogis(1, shape = LL, scale = b);
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
	
	if(SystemType=='<Gamma(G)_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		s = rgamma(N,shape=G,scale = a1/G);
		b = (b1*LL*sin(pi/LL))/pi;
		
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
						if(i==1){sr=rllogis(1, shape = LL, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rgamma(1,shape=G,scale = a1/G);
						snv = c(s_r[s_r>0], s_new);
						sr=rllogis(1, shape = LL, scale = b);
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
	
	if(SystemType=='<Lognormal(sig)_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = log(a1)-(sig^2)/2;
		s = rlnorm(N, meanlog =a, sdlog = sig);
		b = (b1*LL*sin(pi/LL))/pi;
		
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
						if(i==1){sr=rllogis(1, shape = LL, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rlnorm(1, meanlog =a, sdlog = sig);
						snv = c(s_r[s_r>0], s_new);
						sr = rllogis(1, shape = LL, scale = b);
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
	if(SystemType=='<Loglogistic(LL)_k/Expo/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
		
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						s_new = rllogis(1, shape = LL, scale = a);
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

	if(SystemType=='<Loglogistic(LL)_k/Weibull(W)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
		
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						s_new = rllogis(1, shape = LL, scale = a);
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
	
	if(SystemType=='<Loglogistic(LL)_k/Pareto(P)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
		
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						s_new = rllogis(1, shape = LL, scale = a);
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
	
	if(SystemType=='<Loglogistic(LL)_k/Gamma(G)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
		
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						s_new = rllogis(1, shape = LL, scale = a);
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
	
	if(SystemType=='<Loglogistic(LL)_k/Lognormal(sig)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						s_new = rllogis(1, shape = LL, scale = a);
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
	
	if(SystemType=='<Loglogistic(LL)_k/Loglogistic(LL)/1>'){
		s=vector('numeric',N);
		snv=vector('numeric',N);
		tv=vector();
		sr=vector();
		
		a = (a1*LL*sin(pi/LL))/pi;
		s = rllogis(N, shape = LL, scale = a);
		b = (b1*LL*sin(pi/LL))/pi;
		
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
					s_new = rllogis(1, shape = LL, scale = a);
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
						if(i==1){sr=rllogis(1, shape = LL_b, scale = b);}
						else if(length(tv)>2){
							sr = abs(sr[] - tp); 
						}
					}
					
					if(t==t_rem){
						s_new = rllogis(1, shape = LL, scale = a);
						snv = c(s_r[s_r>0], s_new);
						sr=rllogis(1, shape = LL_b, scale = b);
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

	if(SystemType=='<Expo_k/Expo/1>'){	
		num = seq(1,a1);
		val = length(num);
		a1 = vector('numeric',val);
		a = vector('numeric',val);
		L = matrix(0,K,val);
		P_EXP = matrix(0,K,val);
		A = matrix(0,K,val);
		ro = vector('numeric',val);
		SA31 = vector('numeric',val);
		SA32 = vector('numeric',val);
		SA33 = vector('numeric',val);
		SA41 = vector('numeric',val);
		SA42 = vector('numeric',val);
		SA43 = vector('numeric',val);
		SI41 = vector('numeric',val);
		SI42 = vector('numeric',val);
		C1 = matrix(0,K+1,val);
		CI = vector('numeric',val);
		CK = vector('numeric',val);
		C = vector('numeric',val);
		Pk_EXP = vector('numeric',val);
		Pk = vector('numeric',val);
		
		for(j in 1:val){
			a1[j] = j; 
			a[j] = 1/a1[j];
			b = 1/b1;
			ro[j] = a1[j]/b1;
			
			F1 = function(x){dexp(x, rate = b, log = FALSE)} 
			#F1 = function(x){b*exp(-b*x)} 
			#b1 = integrate(function(x){x*F1(x)},0,Inf)$value
			#b1 = 1/b
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_EXP[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},0,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_EXP[2,j])*(1/P_EXP[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_EXP[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_EXP[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_EXP[2,j])*(1/P_EXP[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_EXP[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_EXP[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_EXP[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_EXP[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_EXP[2,j]/L[1,j];
			C1[2,j] = (1-P_EXP[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_EXP[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_EXP[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_EXP[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_EXP[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_EXP[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1
			
			if(K==2){
				Pk_EXP[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_EXP[K,j])/L[K,j]));
			}
			else{
				Pk_EXP[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_EXP[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_EXP[j];
		}
	}

	if(SystemType=='<Expo_k/Weibull(W)/1>'){	
		num = seq(1,a1);
		val = length(num);
		a1 = vector('numeric',val);
		a = vector('numeric',val);
		L = matrix(0,K,val);
		P_WB = matrix(0,K,val);
		A = matrix(0,K,val);
		ro = vector('numeric',val);
		SA31 = vector('numeric',val);
		SA32 = vector('numeric',val);
		SA33 = vector('numeric',val);
		SA41 = vector('numeric',val);
		SA42 = vector('numeric',val);
		SA43 = vector('numeric',val);
		SI41 = vector('numeric',val);
		SI42 = vector('numeric',val);
		C1 = matrix(0,K+1,val);
		CI = vector('numeric',val);
		CK = vector('numeric',val);
		C = vector('numeric',val);
		Pk_WB = vector('numeric',val);
		Pk = vector('numeric',val);
	
		for(j in 1:val){
			a1[j] = j; 
			a[j] = 1/a1[j];
			b = b1/gamma(1+(1/W));
			ro[j] = a1[j]/b1;
			
			F1 = function(x){dweibull(x, shape = W, scale = b, log = FALSE)} 
			#F1 = function(x){(W/b)*((x/b)^(W-1))*exp(-(x/b)^W)}
			#b1 = integrate(function(x){x*F1(x)},0,Inf)$value
			#b1 = b*gamma(1+(1/W));
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_WB[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},0,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_WB[2,j])*(1/P_WB[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_WB[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_WB[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_WB[2,j])*(1/P_WB[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_WB[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_WB[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_WB[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_WB[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_WB[2,j]/L[1,j];
			C1[2,j] = (1-P_WB[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_WB[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_WB[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_WB[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_WB[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_WB[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1;
			
			if(K==2){
				Pk_WB[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_WB[K,j])/L[K,j]));
			}
			else{
				Pk_WB[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_WB[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_WB[j];
		}
	}

	if(SystemType=='<Expo_k/Pareto(P)/1>'){
		
	    num = seq(1,a1);
	    val = length(num);
	    a1 = vector('numeric',val);
	    a = vector('numeric',val);
	    L = matrix(0,K,val);
	    P_PAR = matrix(0,K,val);
	    A = matrix(0,K,val);
	    ro = vector('numeric',val);
	    SA31 = vector('numeric',val);
	    SA32 = vector('numeric',val);
	    SA33 = vector('numeric',val);
	    SA41 = vector('numeric',val);
	    SA42 = vector('numeric',val);
	    SA43 = vector('numeric',val);
	    SI41 = vector('numeric',val);
	    SI42 = vector('numeric',val);
	    C1 = matrix(0,K+1,val);
	    CI = vector('numeric',val);
	    CK = vector('numeric',val);
	    C = vector('numeric',val);
	    Pk_PAR = vector('numeric',val);
	    Pk = vector('numeric',val);
		
	    for(j in 1:val){
			a1[j] = j; 
		    a[j] = 1/a1[j];
			b = b1*(P-1)/P;
			ro[j] = a1[j]/b1;
			
			F1 = function(x){dpareto(x, shape = P, scale = b, log = FALSE)} 
			#F1 = function(x){(P*(b^P))/(x^(P+1))}
			#b1 = integrate(function(x){x*F1(x)},b,Inf)$value
			#b1 = (P*b)/(P-1); #### P>1; else "Inf"
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_PAR[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},b,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_PAR[2,j])*(1/P_PAR[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_PAR[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_PAR[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_PAR[2,j])*(1/P_PAR[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_PAR[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_PAR[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_PAR[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_PAR[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_PAR[2,j]/L[1,j];
			C1[2,j] = (1-P_PAR[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_PAR[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_PAR[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_PAR[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_PAR[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_PAR[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1;
			
			if(K==2){
				Pk_PAR[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_PAR[K,j])/L[K,j]));
			}
			else{
				Pk_PAR[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_PAR[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_PAR[j];
		}
	}

	if(SystemType=='<Expo_k/Gamma(G)/1>'){	
		num = seq(1,a1);
		val = length(num);
		a1 = vector('numeric',val);
		a = vector('numeric',val);
		L = matrix(0,K,val);
		P_G = matrix(0,K,val);
		A = matrix(0,K,val);
		ro = vector('numeric',val);
		SA31 = vector('numeric',val);
		SA32 = vector('numeric',val);
		SA33 = vector('numeric',val);
		SA41 = vector('numeric',val);
		SA42 = vector('numeric',val);
		SA43 = vector('numeric',val);
		SI41 = vector('numeric',val);
		SI42 = vector('numeric',val);
		C1 = matrix(0,K+1,val);
		CI = vector('numeric',val);
		CK = vector('numeric',val);
		C = vector('numeric',val);
		Pk_G = vector('numeric',val);
		Pk = vector('numeric',val);
		
		for(j in 1:val){
			a1[j] = j; 
			a[j] = 1/a1[j];
			b = b1/G;
			ro[j] = a1[j]/b1;
			
			F1 = function(x){dgamma(x, shape = G, scale = b, log = FALSE)} 
			#F1 = function(x){(1/(gamma(G)*b^G))*(x^(G-1) * exp(-x/b))}
			#b1 = integrate(function(x){x*F1(x)},0,Inf)$value
			#b1 = b*G
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_G[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},0,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_G[2,j])*(1/P_G[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_G[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_G[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_G[2,j])*(1/P_G[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_G[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_G[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_G[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_G[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_G[2,j]/L[1,j];
			C1[2,j] = (1-P_G[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_G[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_G[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_G[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_G[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_G[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1;
			
			if(K==2){
				Pk_G[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_G[K,j])/L[K,j]));
			}
			else{
				Pk_G[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_G[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_G[j];
		}
	}

	if(SystemType=='<Expo_k/Lognormal(sig)/1>'){
		num = seq(1,a1);
		val = length(num);
		a1 = vector('numeric',val);
		a = vector('numeric',val);
		L = matrix(0,K,val);
		P_LN = matrix(0,K,val);
		A = matrix(0,K,val);
		ro = vector('numeric',val);
		SA31 = vector('numeric',val);
		SA32 = vector('numeric',val);
		SA33 = vector('numeric',val);
		SA41 = vector('numeric',val);
		SA42 = vector('numeric',val);
		SA43 = vector('numeric',val);
		SI41 = vector('numeric',val);
		SI42 = vector('numeric',val);
		C1 = matrix(0,K+1,val);
		CI = vector('numeric',val);
		CK = vector('numeric',val);
		C = vector('numeric',val);
		Pk_LN = vector('numeric',val);
		
		Pk = vector('numeric',val);
		for(j in 1:val){
			a1[j] = j; 
			a[j] = 1/a1[j];
			b = log(b1)-(sig^2)/2;
			ro[j] = a1[j]/b1;
			
			F1 = function(x){dlnorm(x, meanlog =b, sdlog = sig, log = FALSE)} 
			#F1 = function(x){(1/(x*sig*sqrt(2*pi)))*exp(-((log(x)-b)^2 / (2*sig^2)))}
			#b1 = integrate(function(x){x*F1(x)},0,Inf)$value
			#b1 = exp(b+(sig^2/2))
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_LN[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},0,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_LN[2,j])*(1/P_LN[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_LN[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_LN[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_LN[2,j])*(1/P_LN[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_LN[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_LN[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_LN[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_LN[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_LN[2,j]/L[1,j];
			C1[2,j] = (1-P_LN[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_LN[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_LN[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_LN[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_LN[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_LN[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1;
			
			if(K==2){
				Pk_LN[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_LN[K,j])/L[K,j]));
			}
			else{
				Pk_LN[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_LN[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_LN[j];
		}
	}
	
	if(SystemType=='<Expo_k/Loglogistic(LL)/1>'){	
		num = seq(1,a1);
		val = length(num);
		a1 = vector('numeric',val);
		a = vector('numeric',val);
		L = matrix(0,K,val);
		P_LL = matrix(0,K,val);
		A = matrix(0,K,val);
		ro = vector('numeric',val);
		SA31 = vector('numeric',val);
		SA32 = vector('numeric',val);
		SA33 = vector('numeric',val);
		SA41 = vector('numeric',val);
		SA42 = vector('numeric',val);
		SA43 = vector('numeric',val);
		SI41 = vector('numeric',val);
		SI42 = vector('numeric',val);
		C1 = matrix(0,K+1,val);
		CI = vector('numeric',val);
		CK = vector('numeric',val);
		C = vector('numeric',val);
		Pk_LL = vector('numeric',val);
		
		Pk = vector('numeric',val);
		for(j in 1:val){
			a1[j] = j; 
			a[j] = 1/a1[j];
			b = (b1*LL*sin(pi/LL))/pi;
			ro[j] = a1[j]/b1;
			
            F1 = function(x){dllogis(x, shape = LL, scale = b, log = FALSE)}
            #F1 = function(x){((LL/b)*(x/b)^(LL-1)) / (1 + (x/b)^LL)^2}
			
			for(i in 0:(K-1)){
				L[i+1,j]=(N-i)*a[j];
				P_LL[i+1,j] = integrate(function(x){exp(-L[i+1,j]*x)*F1(x)},0,Inf)$value
			}
			
			if(K==2){
				for(i in 1:K){
					A[1,j] = 1;			
					A[K,j] = A[K-1,j];
				}
			}
			
			if(K==3){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_LL[2,j])*(1/P_LL[3,j]);
					
					SA31[j] = A[K-1,j]*(1 + P_LL[K,j]);
					SA32[j] = 0;
					SA33[j] = 0;
					for(v in 1:(K-2)){
						SA32[j] = SA32[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA33[j] = SA33[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_LL[v+1,j]);
					}
					A[K,j] = SA31[j] + SA32[j] - SA33[j];
				}
			}
			
			if(K>=4){
				for(i in 1:K){
					A[1,j] = 1;			
					A[2,j] = (1-(1-(L[2,j]/(L[2,j]-L[3,j])))*P_LL[2,j])*(1/P_LL[3,j]);
					
					for(v in 2:(K-2)){
						SI41[j] = 0;
						for(w in 1:(v-1)){
							SI41[j] = SI41[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]);
						}
						SI42[j] = 0;
						for(w in 1:v){
							SI42[j] = SI42[j] + sum(((-1)^(v+1-w))*(prod(L[(w:v)+1,j]/(L[w+1,j]-L[(w:v)+2,j])))*A[w,j]*P_LL[w+1,j]);
						}
						A[v+1,j] = (A[v,j] + SI41[j] - SI42[j])*(1/P_LL[v+2,j]);
					}
					
					SA41[j] = A[K-1,j]*(1 + P_LL[K,j]);
					SA42[j] = 0;
					SA43[j] = 0;
					for(v in 1:(K-2)){
						SA42[j] = SA42[j] + sum(((-1)^(K-1-v))*(prod(L[(v:(K-2))+1,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]);
						SA43[j] = SA43[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*P_LL[v+1,j]);
					}
					
					A[K,j] = SA41[j] + SA42[j] - SA43[j];
				}
			}
			
			C1[1,j] = P_LL[2,j]/L[1,j];
			C1[2,j] = (1-P_LL[2,j])/L[2,j];
			
			if(K>2){
				for(v in 2:(K-1)){
					CI[j] = 0;
					for(w in 1:(v-1)){
						CI[j] = CI[j] + sum(((-1)^(v-w))*(prod(L[(w:(v-1))+1,j]/(L[w+1,j]-L[(w:(v-1))+2,j])))*A[w,j]*((1-P_LL[w+1,j])/L[w+1,j]));
					}
					C1[v+1,j] = A[v,j]*((1-P_LL[v+1,j])/L[v+1,j]) + CI[j];
				}
			}
			CK[j] = 0;
			if(K==2){
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_LL[K,j])/L[K,j]);
			}
			else{
				for(v in 1:(K-2)){
					CK[j] = CK[j] + sum(((-1)^(K-v))*(prod(L[(v:(K-2))+2,j]/(L[v+1,j]-L[(v:(K-2))+2,j])))*A[v,j]*((1-P_LL[v+1,j])/L[v+1,j]));
				}
				C1[K+1,j] = A[K,j]*b1 - A[K-1,j]*((1-P_LL[K,j])/L[K,j]) + CK[j];
			}
			
			C[j] = (sum(C1[,j]))^-1;
			
			if(K==2){
				Pk_LL[j] = C[j]*(A[K,j]*b1 - A[K-1,j]*((1-P_LL[K,j])/L[K,j]));
			}
			else{
				Pk_LL[j] = C[j]*((A[K,j]*b1 - A[K-1,j]*((1-P_LL[K,j])/L[K,j]) + CK[j]));
			}
			Pk[j] = Pk_LL[j];
		}
	}

return(Pk) # Data of the graph of the probability of failure-free system (Данные Графика вероятности безотказной работы системы)
}


######################################################################################
######################################################################################

############install.packages("VGAM",dependencies=T,repos='http:// cran.rstudio.com/')
#### install.packages("installr"); require(installr); updateR();

main=function(SystemType=NULL){

	if(is.null(SystemType)==TRUE){
		SystemType=select.list(c('<Expo_k/Expo/1>','<Expo_k/Weibull(W)/1>','<Expo_k/Pareto(P)/1>','<Expo_k/Gamma(G)/1>','<Expo_k/Lognormal(sig)/1>','<Expo_k/Loglogistic(LL)/1>',
		'<Weibull(W)_k/Expo/1>','<Weibull(W)_k/Weibull(W)/1>','<Weibull(W)_k/Pareto(P)/1>','<Weibull(W)_k/Gamma(G)/1>','<Weibull(W)_k/Lognormal(sig)/1>','<Weibull(W)_k/Loglogistic(LL)/1>',
		'<Pareto(P)_k/Expo/1>','<Pareto(P)_k/Weibull(W)/1>','<Pareto(P)_k/Pareto(P)/1>','<Pareto(P)_k/Gamma(G)/1>','<Pareto(P)_k/Lognormal(sig)/1>','<Pareto(P)_k/Loglogistic(LL)/1>',
		'<Gamma(G)_k/Expo/1>','<Gamma(G)_k/Weibull(W)/1>','<Gamma(G)_k/Pareto(P)/1>','<Gamma(G)_k/Gamma(G)/1>','<Gamma(G)_k/Lognormal(sig)/1>','<Gamma(G)_k/Loglogistic(LL)/1>',
		'<Lognormal(sig)_k/Expo/1>','<Lognormal(sig)_k/Weibull(W)/1>','<Lognormal(sig)_k/Pareto(P)/1>','<Lognormal(sig)_k/Gamma(G)/1>','<Lognormal(sig)_k/Lognormal(sig)/1>','<Lognormal(sig)_k/Loglogistic(LL)/1>',
		'<Loglogistic(LL)_k/Expo/1>','<Loglogistic(LL)_k/Weibull(W)/1>','<Loglogistic(LL)_k/Pareto(P)/1>','<Loglogistic(LL)_k/Gamma(G)/1>','<Loglogistic(LL)_k/Lognormal(sig)/1>','<Loglogistic(LL)_k/Loglogistic(LL)/1>'),
		preselect = NULL, multiple = FALSE, title = NULL, graphics = getOption("menu.graphics"))
	}
	
	if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Weibull(W)_k/Expo/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Pareto(P)_k/Expo/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Gamma(G)_k/Expo/1>'||
		SystemType=='<Expo_k/Lognormal(sig)/1>'||SystemType=='<Lognormal(sig)_k/Expo/1>'||SystemType=='<Expo_k/Loglogistic(LL)/1>'||SystemType=='<Loglogistic(LL)_k/Expo/1>'||
		SystemType=='<Weibull(W)_k/Weibull(W)/1>'||SystemType=='<Weibull(W)_k/Pareto(P)/1>'||SystemType=='<Weibull(W)_k/Gamma(G)/1>'||SystemType=='<Weibull(W)_k/Lognormal(sig)/1>'||SystemType=='<Weibull(W)_k/Loglogistic(LL)/1>'||
		SystemType=='<Pareto(P)_k/Weibull(W)/1>'||SystemType=='<Pareto(P)_k/Pareto(P)/1>'||SystemType=='<Pareto(P)_k/Gamma(G)/1>'||SystemType=='<Pareto(P)_k/Lognormal(sig)/1>'||SystemType=='<Pareto(P)_k/Loglogistic(LL)/1>'||
		SystemType=='<Gamma(G)_k/Weibull(W)/1>'||SystemType=='<Gamma(G)_k/Pareto(P)/1>'||SystemType=='<Gamma(G)_k/Gamma(G)/1>'||SystemType=='<Gamma(G)_k/Lognormal(sig)/1>'||SystemType=='<Gamma(G)_k/Loglogistic(LL)/1>'||
		SystemType=='<Lognormal(sig)_k/Weibull(W)/1>'||SystemType=='<Lognormal(sig)_k/Pareto(P)/1>'||SystemType=='<Lognormal(sig)_k/Gamma(G)/1>'||SystemType=='<Lognormal(sig)_k/Lognormal(sig)/1>'||SystemType=='<Lognormal(sig)_k/Loglogistic(LL)/1>'||
		SystemType=='<Loglogistic(LL)_k/Weibull(W)/1>'||SystemType=='<Loglogistic(LL)_k/Pareto(P)/1>'||SystemType=='<Loglogistic(LL)_k/Gamma(G)/1>'||SystemType=='<Loglogistic(LL)_k/Lognormal(sig)/1>'||SystemType=='<Loglogistic(LL)_k/Loglogistic(LL)/1>'){
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
	if(SystemType=='<Expo_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL<<-scan(nlines=1);
	}
	if(SystemType=='<Loglogistic(LL)_k/Expo/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1);
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
	if(SystemType=='<Weibull(W)_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) W "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');W<<-scan(nlines=1);
		cat('input parameter (вводите параметр) LL "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL<<-scan(nlines=1);
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
	if(SystemType=='<Pareto(P)_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) P "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');P<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) LL "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL<<-scan(nlines=1);
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
	if(SystemType=='<Gamma(G)_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) G "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');G<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) LL "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL<<-scan(nlines=1);
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
	if(SystemType=='<Lognormal(sig)_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) sig "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');sig<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) LL "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL<<-scan(nlines=1);
	}
#
	if(SystemType=='<Loglogistic(LL)_k/Weibull(W)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) W "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');W<<-scan(nlines=1);
	}
	if(SystemType=='<Loglogistic(LL)_k/Pareto(P)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) P "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');P<<-scan(nlines=1);
	}
	if(SystemType=='<Loglogistic(LL)_k/Gamma(G)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) G "repair time distribution parameter (параметр распределения времени ремонта)": \n');G<<-scan(nlines=1);
	}
	if(SystemType=='<Loglogistic(LL)_k/Lognormal(sig)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) sig "Parameter of repair time distribution element  (параметр распределения времени ремонта элемента)": \n');sig<<-scan(nlines=1);
	}
	if(SystemType=='<Loglogistic(LL)_k/Loglogistic(LL)/1>'){
		cat('input parameter (вводите параметр) LL "Parameter of failure-free distribution element (параметр распределения времени безотказной работы системы элемента)": \n');LL<<-scan(nlines=1); 
		cat('input parameter (вводите параметр) LL_b "Parameter of repair time distribution element (параметр распределения времени ремонта элемента)": \n');LL_b<<-scan(nlines=1);
	}
#
	if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){ #||SystemType=='<Expo_k/Loglogistic(LL)/1>'){
		Pk_a = SysA_process(SystemType);
	}

###########################################################################
## Trajectory graphics (Графики траекторий) 

	pr=Sys_process(a1, b1, K, N, TT,NG, SystemType);

	NumF=pr[1,3,]
	plot(c(0,pr[1,1,]),c(0,NumF),col="green",type="s", yaxt='n', xaxt='n',xlab="t", ylab="v")
	mtext("t",side=1,at=T+5,line=1,cex=1.1)
	axis(side = 2, at = seq(0,K), las=1)
	axis(side = 1, at = c(0,pr[1,1,]), las=2)
	abline(v=TT,col="red")

################################################
## graph of the probability of failure-free system (График вероятности безотказной работы системы)

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
	pi_K=numeric();
	l=0;
	for(a1 in seq(1,a1)){l=l+1
		for(j in 1:NG){
			
			prv=Sys_process(a1, b1, K, N, TT, NG, SystemType);
			
			pivv=vector('numeric',K+1)
			for(v in 0:(K)){
				pivv[v+1] = summ_T(prv,v)/TT;
			}
			
			viborka[j]=pivv[K+1];
		}
		pi_K[l] = mean(viborka); 
	}

	x11()
	ro = seq(1,a1);
	if(SystemType=='<Expo_k/Expo/1>'||SystemType=='<Expo_k/Weibull(W)/1>'||SystemType=='<Expo_k/Pareto(P)/1>'||SystemType=='<Expo_k/Gamma(G)/1>'||SystemType=='<Expo_k/Lognormal(sig)/1>'){ #||SystemType=='<Expo_k/Loglogistic(LL)/1>'){
		plot(ro,1-pi_K,col="green",type="l",lty=1,lwd=2,xlab=expression(paste(rho," = ",EA/EB)), ylab=expression(1-pi[k]),xaxp = c(1,50,49))
		grid(NULL,lty="dotted")
		legend("right",bg="white",  legend=c('<GI_k|GI|1> Simul.(Имит.)','<Expo_k|GI|1> Analit.(Аналит.)'), lty=c(1,5), col = c("green","red"), inset = .02)
		#legend("right",bg="white",  legend=c(SystemType), inset = .02)
		lines(ro,1-Pk_a,ylim=c(0:1),type='l',ylim=c(0:1),lty=5,col="red",lwd=2)
		
		cat('Vector of system uptime probability by analytical  (Вектор вероятности безотказной работы системы по аналитическому):\n')
		print(1-Pk_a)
		
		cat('Allowable error (Допустимая погрешность):\n')
		PP = abs(pi_K-Pk_a)
		ss = max(PP);ss  
		print(ss)
		
		PFFSO = c(1-Pk_a, 1-pi_K);
	}
	else{
		plot(ro,1-pi_K,col="green",type="l",lwd=2,ylim=c(0:1),lty=1,xlab=expression(paste(rho," = ",EA/EB)), ylab=expression(1-pi[n]),xaxp = c(1,50,49))
		grid(NULL,lty="dotted")
		legend("right",bg="white",  legend=c(SystemType), inset = .02)

		PFFSO = c(1-pi_K);
	}
	
	cat('Vector of system uptime probability by simulation  (Вектор вероятности безотказной работы системы по имитационному):\n')
	return(c(PFFSO)) #Data of the graph of the probability of failure-free system (Данные Графика вероятности безотказной работы системы)
}
#####################################################
#####################################################

#Pk = main();Pk



