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

Sys_process=function(a1, b1, K, N, NG,SystemType=NULL){

	#a1	- Average time between element failures (Среднее время между отказами элементов)
	#b1	- Average repair time (Среднее время ремонта)
	#K	- The number of elements operating in the system (Число элементов обеспечивающих оперативное функционирование системы)
	#N	- The number of elements in the system (Число элементов в системе)
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		b=log(b1)-(sig_b^2)/2;
		s=rlnorm(N, meanlog =a, sdlog = sig);
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
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
		
		while(t<Inf){
			
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
				break;				
			}
			
			r=array(r, dim=c(1,3,k));
			r[,,k] <-c(t,i,j); 
			i=j;
			k=k+1;
		}
	}

	
	return(r)
}

######################################################################################
######################################################################################

############install.packages("VGAM",dependencies=T,repos='hInfp:// cran.rstudio.com/')
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

###########################################################################
##Trajectory Graphs(Графики траекторий)

	pr=Sys_process(a1, b1, K, N, NG,SystemType);

		NumF=pr[1,3,]
		plot(c(0,pr[1,1,]),c(0,NumF),col="green",type="s", lwd=2, yaxt='n', xaxt='n',xlab="t", ylab="v")
		mtext("t",side=1,at=T+5,line=1,cex=1.1)
		axis(side = 2, at = seq(0,K), las=1)
		axis(side = 1, at = c(0,pr[1,1,]), las=2)
		abline(v=max(c(0,pr[1,1,])),col="red", lty=2)

################################################
##Estimation of the Average Lifetime of a System (Оценка среднего времени жизни системы)
	SumT = vector('numeric',NG)
	for(i in 1:NG){

		pr=Sys_process(a1, b1, K, N, NG,SystemType);

		SumT[i] = pr[1,1,dim(pr)[3]];
	}
	#cat('Vector of lifetimes of the system containing NG values (Вектор времен жизни системы содержащий NG значений):\n')
	ET = SumT
	#print(ET)

	MeanR = mean(SumT)

	EPF = ecdf(ET)
	cat('Summary (резюме):\n')
	summary(EPF)
	print(summary(EPF))

	#x = unique(sort(c(seq(0, max(ET), length = 500), knots(EPF))))
	x = unique(sort(c(seq(0, max(ET)))))

	Fx = EPF(x)
	Rx = 1 - Fx

	x11()
	plot(Fx, type="l", col = "blue", main = "Empirical function F*(t) and reliability function R*(t)(Эмпирическая функция F*(t) и функция надёжности R*(t))" , xlab=expression(paste(t)),ylab=expression("F*(t);R*(t)"))
	legend("right",bg="white",  legend=c('F*(t)','R*(t)'), lty=c(1,1), col = c("blue","red"), inset = .02)
	#legend(locator(1),bg="white",  legend=c('F*(t)','R*(t)'), lty=c(1,1), col = c("blue","red"), inset = .02)

	lines(Rx, col = "red")

	cat('Estimation of the average system lifetime (Оценка среднего времени жизни системы):\n')
	return(MeanR)
}

###############################################

#ET=main(); ET