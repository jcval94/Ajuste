ajusteCl<-function(X,n.obs=length(X),ref="OP",crt=1,plot=F,subplot=F,p.val_adms=.05){
  if(!is.numeric(ref)){}
  else{ if(ref>length(X)/3){warning("Number of clusters must be less than input length/3")
    return(NULL)}}
  desc<-function(X,fns=F,ref.=ref,crt.=crt,subplot.=subplot,p.val_adms.=p.val_adms){
    eval<-function(X,fns.=fns,crt.=crt,subplot.=subplot,p.val_adms.=p.val_adms){
      Ajuste<-ajustar1(X,length(X),criteria = crt,plot = subplot,p.val_adms=p.val_adms)
      Ajuste
    }
    div<-function(X,ref.=ref){
      df<-data.frame(A=1:length(X),B=X)
      Enteros<-X-floor(X)==0
      if(any(Enteros)){
        df$CL<-ifelse(Enteros,1,2)
      }else{
        if(!is.numeric(ref)){
          library(mclust)
          mod1<-Mclust(X)$classification
          if(length(table(mod1))==1){
            df$CL<-kmeans(df,2)$cluster
          }else{
            df$CL<-mod1
          }
        }else{
          df$CL<-kmeans(df,ref)$cluster
        }
      }
      CLS<-purrr::map(unique(df$CL),~df[df$CL==.x,2])
      CLS
      return(CLS)
    }
    suppressWarnings(EV<-eval(X,fns))
    if(is.null(EV)){
      if(length(X)>40){
        DV<-purrr::map(div(X),~desc(.x,fns))
        return(DV)
      }else{
        FN<-rnorm
        formals(FN)[1]<-length(X)
        formals(FN)[2]<-mean(X)
        formals(FN)[3]<-ifelse(length(X)==1,0,sd(X))
        return(list(paste("normal(",mean(X),",",ifelse(length(X)==1,0,sd(X)),")"),FN,FN(),data.frame(AD_p.v=0,KS_p.v=0,Chs_p.v=0)))
      }
    }else{
      return(EV)
    }
  }
  FCNS<-desc(X)
  flattenlist <- function(x){  
    morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
    out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
    if(sum(morelists)){ 
      base::Recall(out)
    }else{
      return(out)
    }
  }
  superficie<-flattenlist(FCNS)
  FUN<-superficie[purrr::map_lgl(superficie,~"function" %in% class(.x))]
  Global_FUN<-superficie[purrr::map_lgl(superficie,~"gl_fun" %in% class(.x))]
  Dist<-unlist(superficie[purrr::map_lgl(superficie,is.character)])
  PLTS<-superficie[purrr::map_lgl(superficie,is.ggplot)]
  PV<-do.call("rbind",superficie[purrr::map_lgl(superficie,is.data.frame)])
  Len<-MA<-c()
  repp<-floor(n.obs/length(X))+1
  for (OBS in 1:repp) {
    for (mst in 1:length(FUN)) {
      ljsd<-FUN[[mst]]()
      MA<-c(MA,ljsd)
      if(OBS==1){
        Len<-c(Len,length(ljsd)/length(X))
      }
    }
  }
  MA<-sample(MA,n.obs)
  p.v<-cbind(data.frame(Distribucion=Dist[nchar(Dist)!=0],Prop_dist=Len[nchar(Dist)!=0]),PV)
  library(ggplot2)
  cp<-plt<-c()
  if(plot){
    DF<-rbind(data.frame(A="Ajuste",DT=MA),
              data.frame(A="Real",DT=X))
    plt <- ggplot(DF,aes(x=DT,fill=A)) + geom_density(alpha=0.55)+ggtitle("Dist. Original")
    plt
  }
  TPlts<-c()
  if(subplot){
    library(cowplot)
    cp<-cowplot::plot_grid(plotlist = PLTS, ncol = floor(sqrt(length(PLTS))))
  }
  TPlts<-list(plt,cp)
  return(list(unlist(FUN),MA,p.v,TPlts,Global_FUN))
}