ajustar1<-function(b,generar=1,Cont=T,inputNA,plot=F,
                   p.val_adms=.05,criteria=2,DPQR=T){
  if(missing(inputNA)){b<-na.omit(b)}
  else{b<-ifelse(is.na(b),inputNA,b)}
  if(length(b)==0){
    break()
  }
  if (length(unique(b))<2) {
    fun_g<-function(n=generar){return(rep(b[1],n))}
    return(list(paste0("norm(",b[1],",0)"),fun_g,rep(b[1],generar),data.frame( AD_p.v=1,KS_p.v=1,Chs_p.v=1),NULL))
  }
  if(prod(b==floor(b))==1){
    Cont<-F
  }
  
  if (length(unique(b))==2) {
    p<-length(b[b==unique(b)[1]])/length(b)
    Ber<-function(p.=p,n=generar){
      stats::runif(n) > (1 - p)
    }
    return(list("Ber",Ber,Ber(p)))
  }
  
  DIS<-list(Nombres=c("exp","pois","beta","gamma","lnorm","norm","weibull","nbinom","hyper","cauchy","t"),
            p=c(pexp,ppois,pbeta,pgamma,plnorm,pnorm,pweibull,pnbinom,phyper,pcauchy,pt),
            d=c(dexp,dpois,dbeta,dgamma,dlnorm,dnorm,dweibull,dnbinom,dhyper,dcauchy,dt),
            q=c(qexp,qpois,qbeta,qgamma,qlnorm,qnorm,qweibull,qnbinom,qhyper,qcauchy,qt),
            r=c(rexp,rpois,rbeta,rgamma,rlnorm,rnorm,rweibull,rnbinom,rhyper,rcauchy,rt),
            d_c=c(1,0,1,1,1,1,1,0,0,1,1),
            indicadora=c("0","0","01","0","0","R","0","0","0","R","R")
  )
  DIS<-purrr::map(DIS,~subset(.x, DIS$d_c==as.numeric(Cont)))
  DIS_0<-purrr::map(DIS,~subset(.x, DIS$indicadora=="0"))
  DIS_R<-purrr::map(DIS,~subset(.x, DIS$indicadora=="R"))
  DIS_01<-purrr::map(DIS,~subset(.x, DIS$indicadora=="01"))
  
  if(sum(purrr::map_dbl(DIS_0,~length(.x)))==0){DIS_0<-NULL}
  if(sum(purrr::map_dbl(DIS_R,~length(.x)))==0){DIS_R<-NULL}
  if(sum(purrr::map_dbl(DIS_01,~length(.x)))==0){DIS_01<-NULL}
  bt<-b
  despl<-0
  escala<-1
  eps<-1E-15
  if (sum(b<0)>0){
    if (sum(b<0)/length(b)<0.03){
      bt<-ifelse(b<0,eps,b)
      b_0<-bt
    }
    else{
      b_0<-bt-min(bt)+eps
      despl<- min(bt)
    }
  }else{
    b_0<-bt
  }
  if(max(b)>1){
    escala<-max(bt)
    b_01<-(bt-despl)/(escala-despl)
  }else{
    b_01<-bt
  }
  
  ajustar_b<-function(bt,dist="",Cont.=Cont){
    if(is.null(dist)){break()}
    Disc<-!Cont
    aju<-list()
    if(!dist %in% DIS_01$Nombres){
      suppressWarnings(aju[[1]]<-try(fitdistrplus::fitdist(bt,dist,method = "mle",discrete = Disc),silent = T))
      
    }
    suppressWarnings(aju[[2]]<-try(fitdistrplus::fitdist(bt,dist,method = "mme",discrete = Disc),silent = T))
    suppressWarnings(aju[[3]]<-try(fitdistrplus::fitdist(bt,dist,method = c("mge"),discrete = Disc),silent = T))
    
    suppressWarnings(aju[[4]]<-try(MASS::fitdistr(bt,dist),silent = T))
    if(!assertthat::is.error(aju[[4]])){aju[[4]]$distname<-dist}
    
    if(assertthat::is.error(aju[[1]]) & assertthat::is.error(aju[[2]]) & 
       assertthat::is.error(aju[[3]]) & assertthat::is.error(aju[[4]])){
      return(list())
      break()
    }
    
    funcionales<-!purrr::map_lgl(aju,~assertthat::is.error(.x))
    
    aju<-aju[funcionales]
    
    return(aju)
    
  }
  
  suppressWarnings(try(aju_0<-purrr::map(DIS_0$Nombres,~ajustar_b(b_0,.x)),silent = T))
  suppressWarnings(try(aju_R<-purrr::map(DIS_R$Nombres,~ajustar_b(bt,.x)),silent = T))
  suppressWarnings(try(aju_01<-purrr::map(DIS_01$Nombres,~ajustar_b(b_01,.x)),silent = T))
  
  AAA<-list(aju_0,aju_R,aju_01)
  AAA<-AAA[purrr::map(AAA,~length(.x))!=0]
  
  bts<-list(b_0,bt,b_01)
  num<-0
  Compe<-data.frame()
  
  for (aju_ls in 1:length(AAA)) {
    aju<-AAA[[aju_ls]]
    bs<-bts[[aju_ls]]
    for (comp in 1:length(aju)) {
      for (ress in 1:length(aju[[comp]])) {
        num<-num+1
        if(length(aju[[comp]])!=0){evaluar<-aju[[comp]][[ress]]}
        else{evaluar<-NULL}
        if (is.null(evaluar) | length(evaluar)==0 | 
            c(NA) %in% evaluar$estimate | c(NaN) %in% evaluar$estimate) {next()}
        
        distname<-evaluar$distname
        dist_pfun<-try(get(paste0("p",distname)),silent = T)
        dist_rfun<-try(get(paste0("r",distname)),silent = T)
        
        if(assertthat::is.error(dist_rfun)){next()}
        
        argumentos<-formalArgs(dist_pfun)
        argumentos<-argumentos[argumentos %in% names(evaluar$estimate)]
        num_param<-length(argumentos)
        evaluar$estimate<-evaluar$estimate[names(evaluar$estimate) %in% argumentos]
        
        if(num_param==1){
          EAD<-try(AD<-ADGofTest::ad.test(bs,dist_pfun,evaluar$estimate[1]),silent = T)
          if (Cont) {EKS<-try(KS<-stats::ks.test(bs,dist_pfun,evaluar$estimate[1]),silent = T)}
          else{KS<-data.frame(p.value=0)}
          
          if(assertthat::is.error(EAD) | assertthat::is.error(EKS)){next()}
          if(is.na(EKS$p.value)){next()}
          Chs<-data.frame(p.value=0)
        }
        if(num_param==2){
          
          suppressWarnings(
            Err_pl<-try(AD<-ADGofTest::ad.test(bs,dist_pfun,evaluar$estimate[1],evaluar$estimate[2]),silent = T))
          
          if (assertthat::is.error(Err_pl)) {
            Err_pl<-try(AD<-ADGofTest::ad.test(bs,dist_pfun,evaluar$estimate[1],,evaluar$estimate[2]),silent = T)
          }
          
          if (Cont) {Err_pl2<-try(KS<-stats::ks.test(bs,dist_pfun,evaluar$estimate[1],evaluar$estimate[2]),silent = T)}
          else{KS<-data.frame(p.value=0)}
          
          if(assertthat::is.error(Err_pl) | assertthat::is.error(Err_pl2)){next()}
          if(is.na(Err_pl2$p.value)){next()}
          
          
          suppressWarnings(
            EE_Chs<-try(dst_chsq<-dist_rfun(length(bs),evaluar$estimate[1],evaluar$estimate[2]))
          )
          if(assertthat::is.error(EE_Chs) | prod(is.na(EE_Chs))==1){
            dst_chsq<-dist_rfun(length(bs),evaluar$estimate[1],,evaluar$estimate[2])
          }
          
          #Chs<-chisq.test(bs,dst_chsq,simulate.p.value = 800)
          Chs<-data.frame(p.value=0)
        }
        pvvv<-p.val_adms
        if(criteria==1){
          crit<-AD$p.value>pvvv | KS$p.value>pvvv | Chs$p.value >pvvv
        }
        else{
          crit<-AD$p.value>(pvvv) & KS$p.value>(pvvv)
        }
        if(crit){
          if(aju_ls %in% 3){
            estimate3=despl
            estimate4=escala
          }
          else if(aju_ls==1){
            estimate3=despl
            estimate4=1
          }
          else{
            estimate3=0
            estimate4=1
          }
          
          Compe<-rbind(Compe,data.frame(Dist=distname,AD_p.v=AD$p.value,KS_p.v=KS$p.value,
                                        Chs_p.v=Chs$p.value,
                                        estimate1=evaluar$estimate[1],estimate2=evaluar$estimate[2],
                                        estimateLL1=estimate3,estimateLL2=estimate4
          ))
        }
        else{
          next()
        }
        
      }
    }
  }
  
  
  if (nrow(Compe)==0) {
    warning("No hubo ajuste")
    return(NULL)
    break()
  }
  Compe$PV_S<-rowSums(Compe[,2:4])
  WNR<-Compe[Compe$PV_S %in% max(Compe$PV_S),][1,]
  
  distW<-WNR$Dist
  paramsW<-WNR[1,names(Compe)[startsWith(names(Compe),"estim")]]
  paramsW<-paramsW[,!is.na(paramsW)]
  if(generar<=0){generar<-1}
  
  generadora_r<-function(n=generar,dist=distW,params=paramsW){
    fn<-get(paste0("r",dist))
    formals(fn)[1]<-n
    for (pr in 1:(length(params)-2)) {
      formals(fn)[pr+1]<-as.numeric(params[pr])
    }
    
    fn()*params[,length(params)]+params[,length(params)-1]
  }
  
  if(DPQR){
    generadoras<-function(x,tipo,dist=distW,params=paramsW){
      fn<-get(paste0(tipo,dist))
      formals(fn)[1]<-x
      for (pr in 1:(length(params)-2)) {
        formals(fn)[pr+1]<-as.numeric(params[pr])
      }
      class(fn)<-"gl_fun"
      fn
    }
    rajuste<<-generadora_r
    class(rajuste)<<-"gl_fun"
    pajuste<<-generadoras(1,"p")
    qajuste<<-generadoras(1,"q")
    dajuste<<-generadoras(1,"d")
  }
  
  MA<-generadora_r()
  
  paramsW2<-purrr::map_df(paramsW,~round(.x,3))
  
  if(paramsW2[,length(paramsW2)]!=1 | paramsW2[,length(paramsW2)-1]!=0){
    distribu<-paste0(WNR$Dist,"(",paste0(paramsW2[,1:(length(paramsW2)-2)],collapse = ", "),")*",paramsW2[,length(paramsW2)],"+",paramsW2[,length(paramsW2)-1])
  }
  else{
    distribu<-paste0(WNR$Dist,"(",paste0(paramsW2[,1:(length(paramsW2)-2)],collapse = ", "),")")
  }
  p<-c()
  library(ggplot2)
  if(plot){
    DF<-rbind(data.frame(A="Ajuste",DT=MA),
              data.frame(A="Real",DT=b))
    
    p <- ggplot(DF,aes(x=DT,fill=A)) + geom_density(alpha=0.4) +ggtitle(distribu)
  }
  
  return(list(distribu,generadora_r,MA,WNR[,2:4],p,list(rajuste,pajuste,dajuste,qajuste)))
}