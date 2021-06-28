library(shiny)
library(shinyIncubator)
library(googleVis)

shinyServer(function(input,output,session){
  
  
  
  ### Disclaimer
  
  
  
  output$disclaimer<-renderText(function(){
    return(c('<script type="text/javascript">',"alert('Data is processed in real time. After you quit the website, your data will NOT be stored nor accessible to anyone.');","</script>"))
  })
  
  
  
  ### Fine tuning
  
  
  
  output$adjust<-renderText(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    return(c("<br>","<br>","Some <b>fine tuning</b> (<b>optional</b>, for advanced users only):","<br>","<br>","Feeding types designation may be ambiguous. Default values are based on the classification of Yeates et al. (1993) and subsequently obtained data. If you think differently, use the drop-down lists below to perform adjustment.","<br><br>"))
  })
  
  
  
  ### Compost + Noticed a bug? Mailto!
  
  
  
  output$mailto<-renderText(function(){
    return(c('If you wish to calculate <b>Compost Maturity Index</b> without providing a taxonomic input file, ','<a href="https://sieriebriennikov.shinyapps.io/compost/">click here</a>','.</br></br>',
             "Noticed a bug? This is a beta version, please be so kind to ",'<a href="mailto:sarasm@inia.es">report any issues.</a>'))
  })
  
  
  
  ### [MI family indices]:sigmaMI
  
  
  
  output$sigMI<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating sigma MI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting c-p values from the database
    
    c<-ncol(taxa)
    i<-1
    cpppvector<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      cppp<-dat[dat["taxon"]==taxon][2]
      cpppvector<-append(cpppvector,cppp)
    }
    
    cpppvector<-as.numeric(cpppvector)
    taxacppp<-rbind(taxa,cpppvector)
    
    i<-1
    sigmivector<-vector()
    
    for (i in 1:r)
    {
      sigmi<-sum(taxacppp[i,]*taxacppp[r+1,])/sum(taxacppp[i,])
      sigmivector<-append(sigmivector,sigmi)
      i<-i+1
    }
    
    plot(treats,sigmivector,ylab="Sigma Maturity Index",main="Sigma Maturity index",las=2)
    
  })
  
  
  
  ### [MI family indices]:MI
  
  
  
  output$MI<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating MI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        feed<-substr(feed,1,1)
      }
      feedvector<-append(feedvector,feed)
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    taxacp<-rbind(taxa,cpvector)
    r<-nrow(taxa)
    taxacp<-subset(taxacp,select=which(taxacp[r+1,]>0))
    
    taxami<-taxacp
    
    i<-1
    mivector<-vector()
    
    for (i in 1:r)
    {
      mi<-sum(taxami[i,]*taxami[r+1,])/sum(taxami[i,])
      mivector<-append(mivector,mi)
      i<-i+1
    }
    
    if (length(which(mivector>0))==0) {
      plot(0,type="n",axes=F,xlab="",ylab="")
      text(x=0.5,labels="MI and MI2-5 values could not be plotted as free-living nematodes are absent in your sample.")
      text(x=0,labels="Scroll down to see the PPI plot.")
    } else {
      plot(treats,mivector,ylab="Maturity Index",main="Maturity index",las=2)
    }
    
  })
  
  
  
  ### [MI family indices]:MI2-5
  
  
  
  output$MI25<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating MI 2-5',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        feed<-substr(feed,1,1)
      }
      feedvector<-append(feedvector,feed)
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    taxacp<-rbind(taxa,cpvector)
    r<-nrow(taxa)
    taxacp<-subset(taxacp,select=which(taxacp[r+1,]>0))
    taxami25<-subset(taxacp,select=which(taxacp[r+1,]>1))
    
    i<-1
    mi25vector<-vector()
    
    for (i in 1:r)
    {
      mi25<-sum(taxami25[i,]*taxami25[r+1,])/sum(taxami25[i,])
      mi25vector<-append(mi25vector,mi25)
      i<-i+1
    }
    
    if (length(which(mi25vector>0))==0) {
      return(NULL)
    } else {
      plot(treats,mi25vector,ylab="Maturity Index 2-5", main="Maturity Index 2-5",las=2)
    }
    
  })
  
  
  
  ### [MI family indices]:PPI
  
  
  
  output$PPI<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating PPI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        feed<-substr(feed,1,1)
      }
      feedvector<-append(feedvector,feed)
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    feedvector<-as.numeric(feedvector)
    
    # extracting p-p values from the database
    
    c<-ncol(taxa)
    i<-1
    ppvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        taxon<-names(taxa)[i]
        pp<-dat[dat["taxon"]==taxon][2]
      } else {
        pp<-0
      }
      ppvector<-append(ppvector,pp)
    }
    
    
    ppvector<-as.numeric(ppvector)
    taxapp<-rbind(taxa,ppvector)
    taxapp<-subset(taxapp,select=which(taxapp[r+1,]>0))
    
    r<-nrow(taxa)
    
    taxappi<-taxapp
    
    i<-1
    ppivector<-vector()
    
    for (i in 1:r)
    {
      ppi<-sum(taxappi[i,]*taxappi[r+1,])/sum(taxappi[i,])
      ppivector<-append(ppivector,ppi)
      i<-i+1
    }
    
    if (length(which(ppivector>0))==0) {
      plot(0,type="n",axes=F,xlab="",ylab="")
      text(x=0.5,labels="PPI values could not be plotted as herbivores are absent in your sample.")
      text(x=0,labels="Scroll down to see the c-p triangle.")
    } else {
      plot(treats,ppivector,ylab="Plant Parasitic Index",main="Plant Parasitic Index",las=2)
    }
    
    
  })
  
  
  
  ### [MI family indices]:c-p triangle
  
  
  
  output$triangle<-renderImage(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Building a c-p triangle',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    outfile<-tempfile(fileext='.png')
    png(outfile,width=900,height=588)
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        feed<-substr(feed,1,1)
      }
      feedvector<-append(feedvector,feed)
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    taxacp<-rbind(taxa,cpvector)
    r<-nrow(taxa)
    taxacp<-subset(taxacp,select=which(taxacp[r+1,]>0))
    
    # drawing the triangle
    
    r<-nrow(taxa)
    
    x<-function(cp1,cp35) 100-cp35-(0.5*cp1)
    y<-function(cp1) 0.5*cp1*sqrt(3)
    
    plot(0,axes=F,xlim=c(-3,150),ylim=c(-10,90),xlab="",ylab="")
    
    points(x(0,0),y(0))
    points(x(100,0),y(100))
    segments(x(0,100),y(0),x(0,0),y(0),col="green",lwd=2.5)
    segments(x(0,0),y(0),x(100,0),y(100),col="blue",lwd=2.5)
    segments(x(100,0),y(100),x(0,100),y(0),col="red",lwd=2.5)
    
    segments(x(0,80),y(0),x(20,80),y(20),col="green")
    segments(x(0,60),y(0),x(40,60),y(40),col="green")
    segments(x(0,40),y(0),x(60,40),y(60),col="green")
    segments(x(0,20),y(0),x(80,20),y(80),col="green")
    
    segments(x(80,0),y(80),x(0,80),y(0),col="blue")
    segments(x(60,0),y(60),x(0,60),y(0),col="blue")
    segments(x(40,0),y(40),x(0,40),y(0),col="blue")
    segments(x(20,0),y(20),x(0,20),y(0),col="blue")
    
    segments(x(20,80),y(20),x(20,0),y(20),col="red")
    segments(x(40,60),y(40),x(40,0),y(40),col="red")
    segments(x(60,40),y(60),x(60,0),y(60),col="red")
    segments(x(80,20),y(80),x(80,0),y(80),col="red")
    
    text(50,-8,labels="c-p 3-5, %",cex=1.6)
    text(8,-8,labels="100 (stability)",cex=1.6)
    text(99,-8,labels="0",cex=1.6)
    
    text(87,44,labels="c-p 2, %",cex=1.6)
    text(113,3,labels="100 (stress)",cex=1.6)
    text(57,90,labels="0",cex=1.6)
    
    text(13,44,labels="c-p 1, %",cex=1.6)
    text(33,90,labels="(enrichment) 100",cex=1.6)
    text(-3,3,labels="0",cex=1.6)
    
    # plotting treatments on the triangle
    
    i<-1
    taxacp1<-subset(taxacp,select=which(taxacp[r+1,]==1))
    taxacp35<-subset(taxacp,select=which(taxacp[r+1,]>2))
    
    symb<-c("A","B","C","D","E","F","G","H","I",'J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','1','2','3','4','5','6','7','8','9','0','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',"AA","BB","CC","DD","EE","FF","GG","HH","II",'JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT','UU','VV','WW','XX','YY','ZZ','11','22','33','44','55','66','77','88','99','00','aa','bb','cc','dd','ee','ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','ww','xx','yy','zz')
    
    for (i in 1:r)
    {
      cp1<-100*sum(taxacp1[i,])/sum(taxacp[i,])
      cp35<-100*sum(taxacp35[i,])/sum(taxacp[i,])
      if (length(levels(treats))<=25) {
        points(x(cp1,cp35),y(cp1),pch=which(levels(treats)==treats[i]),lwd=1.5,cex=1.3)
      } else {
        points(x(cp1,cp35),y(cp1),pch=symb[which(levels(treats)==treats[i])],lwd=1.5,cex=1.3)
      }
      i<-i+1
    }
    
    if (length(levels(treats))<=25) {
      pchseq<-seq(from=1,to=length(levels(treats)))
    } else {
      pchseq<-symb[1:length(levels(treats))]
    }
    
    if (length(levels(treats))<=35) {fontsize<-1} else {fontsize<-0.7}
    legend("topright",legend=levels(treats),pch=pchseq,y.intersp=0.8,cex=fontsize)
    
    dev.off()
    list(src=outfile,
         contentType='image/png',
         width=900,
         height=588,
         alt="C-P triangle")
  },deleteFile=TRUE)
  
  
  
  ### [food web diagnostics]:BI
  
  
  
  output$BI<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating BI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # extracting food web data from the input
    
    feedvectorcopy<-feedvector
    feedvectorcopy[which(feedvector==1)]<-"PF"
    feedvectorcopy[which(feedvector==2)]<-"FF"
    feedvectorcopy[which(feedvector==3)]<-"BF"
    feedvectorcopy[which(feedvector==5)]<-"PRED"
    feedvectorcopy[which(feedvector==7)]<-"PAR"
    feedvectorcopy[which(feedvector==8)]<-"OV"
    groupvector<-paste(feedvectorcopy,cpvector,sep="")
    
    c<-ncol(taxa)
    i<-1
    fwdvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        fwd<-0
      } else {
        if (cpvector[i]==1 || cpvector[i]==4) {
          fwd<-3.2
        }
        if (cpvector[i]==2) {
          fwd<-0.8
        }
        
        if (cpvector[i]==3) {
          fwd<-1.8
        }
        if (cpvector[i]==5) {
          fwd<-5.0
        }
      }
      fwdvector<-append(fwdvector,fwd)
    }
    
    taxafw<-rbind(taxa,groupvector,fwdvector)
    
    r<-nrow(taxa)
    
    taxafw<-subset(taxafw,select=which(taxafw[r+2,]>0))
    
    pred1<-which(taxafw[r+1,]=="PRED1")
    pred2<-which(taxafw[r+1,]=="PRED2")
    pred3<-which(taxafw[r+1,]=="PRED3")
    pred4<-which(taxafw[r+1,]=="PRED4")
    pred5<-which(taxafw[r+1,]=="PRED5")
    ov1<-which(taxafw[r+1,]=="OV1")
    ov2<-which(taxafw[r+1,]=="OV2")
    ov3<-which(taxafw[r+1,]=="OV3")
    ov4<-which(taxafw[r+1,]=="OV4")
    ov5<-which(taxafw[r+1,]=="OV5")
    ff1<-which(taxafw[r+1,]=="FF1")
    ff2<-which(taxafw[r+1,]=="FF2")
    ff3<-which(taxafw[r+1,]=="FF3")
    ff4<-which(taxafw[r+1,]=="FF4")
    ff5<-which(taxafw[r+1,]=="FF5")
    bf1<-which(taxafw[r+1,]=="BF1")
    bf2<-which(taxafw[r+1,]=="BF2")
    bf3<-which(taxafw[r+1,]=="BF3")
    bf4<-which(taxafw[r+1,]=="BF4")
    bf5<-which(taxafw[r+1,]=="BF5")
    pf1<-which(taxafw[r+1,]=="PF1")
    pf2<-which(taxafw[r+1,]=="PF2")
    pf3<-which(taxafw[r+1,]=="PF3")
    pf4<-which(taxafw[r+1,]=="PF4")
    pf5<-which(taxafw[r+1,]=="PF5")
    
    e<-c(bf1,ff2,pred1,ov1)
    taxae<-subset(taxafw,select=e)
    
    b<-c(bf2,ff2)
    taxab<-subset(taxafw,select=b)
    
    s<-c(pred2,pred3,pred4,pred5,ov2,ov3,ov4,ov5,ff3,ff4,ff5,bf3,bf4,bf5)
    taxas<-subset(taxafw,select=s)
    
    i<-1
    bivector<-vector()
    
    for (i in 1:r)
    {
      bi<-100*sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,]))/(sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,]))+sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,])))
      bivector<-append(bivector,bi)
      i<-i+1
    }
    
    if (length(which(bivector>0))==0) {
      plot(0,type="n",axes=F,xlab="",ylab="")
      text(x=0.5,labels="BI values could not be plotted as free-living nematodes are absent in your sample.")
      text(x=0,labels="Switch back to previous tabs.")
    } else {
      plot(treats,bivector,ylab="Basal Index",main="Basal Index",las=2)
    }
    
  })
  
  
  
  ### [food web diagnostics]:CI
  
  
  
  output$CI<-renderPlot(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating CI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # extracting food web data from the input
    
    feedvectorcopy<-feedvector
    feedvectorcopy[which(feedvector==1)]<-"PF"
    feedvectorcopy[which(feedvector==2)]<-"FF"
    feedvectorcopy[which(feedvector==3)]<-"BF"
    feedvectorcopy[which(feedvector==5)]<-"PRED"
    feedvectorcopy[which(feedvector==7)]<-"PAR"
    feedvectorcopy[which(feedvector==8)]<-"OV"
    groupvector<-paste(feedvectorcopy,cpvector,sep="")
    
    c<-ncol(taxa)
    i<-1
    fwdvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        fwd<-0
      } else {
        if (cpvector[i]==1 || cpvector[i]==4) {
          fwd<-3.2
        }
        if (cpvector[i]==2) {
          fwd<-0.8
        }
        
        if (cpvector[i]==3) {
          fwd<-1.8
        }
        if (cpvector[i]==5) {
          fwd<-5.0
        }
      }
      fwdvector<-append(fwdvector,fwd)
    }
    
    taxafw<-rbind(taxa,groupvector,fwdvector)
    
    r<-nrow(taxa)
    
    taxafw<-subset(taxafw,select=which(taxafw[r+2,]>0))
    
    ff2<-which(taxafw[r+1,]=="FF2")
    bf1<-which(taxafw[r+1,]=="BF1")
    
    taxaff2<-subset(taxafw,select=ff2)
    taxabf1<-subset(taxafw,select=bf1)
    
    i<-1
    civector<-vector()
    
    for (i in 1:r)
    {
      ci<-100*sum(as.numeric(taxaff2[r+2,])*as.numeric(taxaff2[i,]))/(sum(as.numeric(taxaff2[r+2,])*as.numeric(taxaff2[i,]))+sum(as.numeric(taxabf1[r+2,])*as.numeric(taxabf1[i,])))
      civector<-append(civector,ci)
      i<-i+1
    }
    
    if (length(which(is.na(civector)==F))==0) {
      plot(0,type="n",axes=F,xlab="",ylab="")
      text(x=0.5,labels="CI values could not be plotted as free-living nematodes are absent in your sample.")
      text(x=0,labels="Switch back to previous tabs.")
    } else {
      plot(treats,civector,ylab="Channel Index",main="Channel Index",las=2)
    }
    
  })
  
  
  
  ### [Food web diagnostics]:FWD square
  
  
  
  output$square<-renderImage(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Conducting food web analysis',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    outfile<-tempfile(fileext=".png")
    png(outfile,width=800,height=650)
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # extracting food web data from the input
    
    feedvectorcopy<-feedvector
    feedvectorcopy[which(feedvector==1)]<-"PF"
    feedvectorcopy[which(feedvector==2)]<-"FF"
    feedvectorcopy[which(feedvector==3)]<-"BF"
    feedvectorcopy[which(feedvector==5)]<-"PRED"
    feedvectorcopy[which(feedvector==7)]<-"PAR"
    feedvectorcopy[which(feedvector==8)]<-"OV"
    groupvector<-paste(feedvectorcopy,cpvector,sep="")
    
    c<-ncol(taxa)
    i<-1
    fwdvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        fwd<-0
      } else {
        if (cpvector[i]==1 || cpvector[i]==4) {
          fwd<-3.2
        }
        if (cpvector[i]==2) {
          fwd<-0.8
        }
        
        if (cpvector[i]==3) {
          fwd<-1.8
        }
        if (cpvector[i]==5) {
          fwd<-5.0
        }
      }
      fwdvector<-append(fwdvector,fwd)
    }
    
    taxafw<-rbind(taxa,groupvector,fwdvector)
    
    r<-nrow(taxa)
    
    taxafw<-subset(taxafw,select=which(taxafw[r+2,]>0))
    
    pred1<-which(taxafw[r+1,]=="PRED1")
    pred2<-which(taxafw[r+1,]=="PRED2")
    pred3<-which(taxafw[r+1,]=="PRED3")
    pred4<-which(taxafw[r+1,]=="PRED4")
    pred5<-which(taxafw[r+1,]=="PRED5")
    ov1<-which(taxafw[r+1,]=="OV1")
    ov2<-which(taxafw[r+1,]=="OV2")
    ov3<-which(taxafw[r+1,]=="OV3")
    ov4<-which(taxafw[r+1,]=="OV4")
    ov5<-which(taxafw[r+1,]=="OV5")
    ff1<-which(taxafw[r+1,]=="FF1")
    ff2<-which(taxafw[r+1,]=="FF2")
    ff3<-which(taxafw[r+1,]=="FF3")
    ff4<-which(taxafw[r+1,]=="FF4")
    ff5<-which(taxafw[r+1,]=="FF5")
    bf1<-which(taxafw[r+1,]=="BF1")
    bf2<-which(taxafw[r+1,]=="BF2")
    bf3<-which(taxafw[r+1,]=="BF3")
    bf4<-which(taxafw[r+1,]=="BF4")
    bf5<-which(taxafw[r+1,]=="BF5")
    
    e<-c(bf1,ff2,pred1,ov1)
    taxae<-subset(taxafw,select=e)
    
    b<-c(bf2,ff2)
    taxab<-subset(taxafw,select=b)
    
    s<-c(pred2,pred3,pred4,pred5,ov2,ov3,ov4,ov5,ff3,ff4,ff5,bf3,bf4,bf5)
    taxas<-subset(taxafw,select=s)
    
    par(type="s",mar=c(5,4,2,12),xpd=TRUE)
    plot(50,50,pch=3,xlim=c(0,100),ylim=c(0,100),xlab="Structure Index",ylab="Enrichment Index",main="Food web analysis",axes=FALSE)
    axis(1,pos=0)
    axis(2,pos=0)
    segments(0,100,100,100)
    segments(100,0,100,100)
    segments(50,0,50,100)
    segments(0,50,100,50)
    
    symb<-c("A","B","C","D","E","F","G","H","I",'J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','1','2','3','4','5','6','7','8','9','0','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',"AA","BB","CC","DD","EE","FF","GG","HH","II",'JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT','UU','VV','WW','XX','YY','ZZ','11','22','33','44','55','66','77','88','99','00','aa','bb','cc','dd','ee','ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','ww','xx','yy','zz')
    
    i<-1
    
    for (i in 1:r)
    {
      si<-100*sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))/(sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      ei<-100*sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))/(sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      if (length(levels(treats))<=25) {
        points(si,ei,pch=which(levels(treats)==treats[i]),lwd=1.5,cex=1.3)
      } else {
        points(si,ei,pch=symb[which(levels(treats)==treats[i])],lwd=1.5,cex=1.3)
      }
      i<-i+1
    }
    
    if (length(levels(treats))<=25) {
      pchseq<-seq(from=1,to=length(levels(treats)))
    } else {
      pchseq<-symb[1:length(levels(treats))]
    }
    if (length(levels(treats))<=35) {fontsize<-1} else {fontsize<-0.7}
    legend("topright",inset=c(-0.15,0),legend=levels(treats),pch=pchseq,y.intersp=0.8,cex=fontsize)
    
    dev.off()
    list(src=outfile,
         contentType='image/png',
         width=800,
         height=650,
         alt="Footprints")
  }, deleteFile=TRUE)
  
  
  
  ### [Metabolic footprints]:Footprint square
  
  
  
  output$footprint<-renderImage(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating metabolic footprints',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    outfile<-tempfile(fileext=".png")
    png(outfile,width=800,height=650)
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # extracting food web data from the input
    
    feedvectorcopy<-feedvector
    feedvectorcopy[which(feedvector==1)]<-"PF"
    feedvectorcopy[which(feedvector==2)]<-"FF"
    feedvectorcopy[which(feedvector==3)]<-"BF"
    feedvectorcopy[which(feedvector==5)]<-"PRED"
    feedvectorcopy[which(feedvector==7)]<-"PAR"
    feedvectorcopy[which(feedvector==8)]<-"OV"
    groupvector<-paste(feedvectorcopy,cpvector,sep="")
    
    c<-ncol(taxa)
    i<-1
    fwdvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        fwd<-0
      } else {
        if (cpvector[i]==1 || cpvector[i]==4) {
          fwd<-3.2
        }
        if (cpvector[i]==2) {
          fwd<-0.8
        }
        
        if (cpvector[i]==3) {
          fwd<-1.8
        }
        if (cpvector[i]==5) {
          fwd<-5.0
        }
      }
      fwdvector<-append(fwdvector,fwd)
    }
    
    i<-1
    footvector<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      mass<-as.numeric(dat[dat["taxon"]==taxon][4])
      cp<-as.numeric(dat[dat["taxon"]==taxon][2])
      foot<-0.1*(mass/cp)+0.273*(mass^0.75)
      footvector<-append(footvector,foot)
      i<-i+1
    }
    
    fwdvector<-as.numeric(fwdvector)
    taxafw<-rbind(taxa,groupvector,fwdvector,footvector)
    
    r<-nrow(taxa)
    
    taxafw<-subset(taxafw,select=which(taxafw[r+2,]>0))
    
    pred1<-which(taxafw[r+1,]=="PRED1")
    pred2<-which(taxafw[r+1,]=="PRED2")
    pred3<-which(taxafw[r+1,]=="PRED3")
    pred4<-which(taxafw[r+1,]=="PRED4")
    pred5<-which(taxafw[r+1,]=="PRED5")
    ov1<-which(taxafw[r+1,]=="OV1")
    ov2<-which(taxafw[r+1,]=="OV2")
    ov3<-which(taxafw[r+1,]=="OV3")
    ov4<-which(taxafw[r+1,]=="OV4")
    ov5<-which(taxafw[r+1,]=="OV5")
    ff1<-which(taxafw[r+1,]=="FF1")
    ff2<-which(taxafw[r+1,]=="FF2")
    ff3<-which(taxafw[r+1,]=="FF3")
    ff4<-which(taxafw[r+1,]=="FF4")
    ff5<-which(taxafw[r+1,]=="FF5")
    bf1<-which(taxafw[r+1,]=="BF1")
    bf2<-which(taxafw[r+1,]=="BF2")
    bf3<-which(taxafw[r+1,]=="BF3")
    bf4<-which(taxafw[r+1,]=="BF4")
    bf5<-which(taxafw[r+1,]=="BF5")
    
    e<-c(bf1,ff2,pred1,ov1)
    taxae<-subset(taxafw,select=e)
    
    b<-c(bf2,ff2)
    taxab<-subset(taxafw,select=b)
    
    s<-c(pred2,pred3,pred4,pred5,ov2,ov3,ov4,ov5,ff3,ff4,ff5,bf3,bf4,bf5)
    taxas<-subset(taxafw,select=s)
    
    i<-1
    footevector<-vector()
    
    for (i in 1:r)
    {
      foote<-sum(as.numeric(taxae[i,])*as.numeric(taxae[r+3,]))
      footevector<-append(footevector,foote)
      i<-i+1
    }
    
    i<-1
    footsvector<-vector()
    
    for (i in 1:r)
    {
      foots<-sum(as.numeric(taxas[i,])*as.numeric(taxas[r+3,]))
      footsvector<-append(footsvector,foots)
      i<-i+1
    }
    
    k<-max(append(footsvector,footevector))/20
    
    i<-1
    sivector<-vector()
    eivector<-vector()
    
    for (i in 1:r)
    {
      si<-100*sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))/(sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      ei<-100*sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))/(sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      sivector<-append(sivector,si)
      eivector<-append(eivector,ei)
      i<-i+1
    }
    
    # plotting the enrichment and structure footprints
    
    matrixformeans<-cbind(eivector,sivector,footevector,footsvector)
    means<-means.by(matrixformeans,treats)
    sds<-sd.by(matrixformeans,treats)
    eimeans<-means[,1]
    simeans<-means[,2]
    footemeans<-means[,3]
    footsmeans<-means[,4]
    footesds<-sds[,3]
    footssds<-sds[,4]
    
    par(mar=c(5,4,2,12),xpd=TRUE)
    plot(50,50,pch=3,xlim=c(0,100),ylim=c(0,100),xlab="Structure Index",ylab="Enrichment Index",main="Metabolic footprint",sub="full line - treatment mean, dotted line - SD",axes=FALSE)
    axis(1,pos=0)
    axis(2,pos=0)
    segments(0,100,100,100)
    segments(100,0,100,100)
    segments(50,0,50,100)
    segments(0,50,100,50)
    
    symb<-c("A","B","C","D","E","F","G","H","I",'J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','1','2','3','4','5','6','7','8','9','0','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z',"AA","BB","CC","DD","EE","FF","GG","HH","II",'JJ','KK','LL','MM','NN','OO','PP','QQ','RR','SS','TT','UU','VV','WW','XX','YY','ZZ','11','22','33','44','55','66','77','88','99','00','aa','bb','cc','dd','ee','ff','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','ww','xx','yy','zz')
    
    i<-1
    l<-length(levels(treats))
    
    for (i in 1:l)
    {
      si<-simeans[i]
      ei<-eimeans[i]
      if (length(levels(treats))<=25) {
        points(si,ei,pch=i,lwd=1.5,cex=1.3)
      } else {
        points(si,ei,pch=symb[i],lwd=1.5)
      }
      
      x1<-si-footsmeans[i]/k
      y1<-ei
      x2<-si
      y2<-ei+footemeans[i]/k
      x3<-si+footsmeans[i]/k
      y3<-ei
      x4<-si
      y4<-ei-footemeans[i]/k
      points(c(x1,x2,x3,x4,x1),c(y1,y2,y3,y4,y1),type="l",lwd=2)
      
      x1<-si-footsmeans[i]/k-footssds[i]/k
      y1<-ei
      x2<-si
      y2<-ei+footemeans[i]/k+footesds[i]/k
      x3<-si+footsmeans[i]/k+footssds[i]/k
      y3<-ei
      x4<-si
      y4<-ei-footemeans[i]/k-footesds[i]/k
      points(c(x1,x2,x3,x4,x1),c(y1,y2,y3,y4,y1),type="l",lty=3)
      
      x1<-si-footsmeans[i]/k+footssds[i]/k
      y1<-ei
      x2<-si
      y2<-ei+footemeans[i]/k-footesds[i]/k
      x3<-si+footsmeans[i]/k-footssds[i]/k
      y3<-ei
      x4<-si
      y4<-ei-footemeans[i]/k+footesds[i]/k
      points(c(x1,x2,x3,x4,x1),c(y1,y2,y3,y4,y1),type="l",lty=3)
      
      i<-i+1
    }
    
    if (length(levels(treats))<=25) {
      pchseq<-seq(from=1,to=length(levels(treats)))
    } else {
      pchseq<-symb[1:length(levels(treats))]
    }
    if (length(levels(treats))<=35) {fontsize<-1} else {fontsize<-0.7}
    legend("topright",inset=c(-0.15,0),legend=levels(treats),pch=pchseq,y.intersp=0.8,cex=fontsize)
    
    dev.off()
    list(src=outfile,
         contentType='image/png',
         width=800,
         height=650,
         alt="Footprints")
  }, deleteFile=TRUE)
  
  
  
  ### [Summary]:Numeric summary table
  
  
  
  output$summary<-renderTable(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculation in progress',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    cpppvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      taxon<-names(taxa)[i]
      cppp<-dat[dat["taxon"]==taxon][2]
      cpppvector<-append(cpppvector,cppp)
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    cpppvector<-as.numeric(cpppvector)
    taxacp<-rbind(taxa,cpvector)
    r<-nrow(taxa)
    taxacp<-subset(taxacp,select=which(taxacp[r+1,]>0))
    taxacppp<-rbind(taxa,cpppvector)
    
    # calculating sigmaMI
    
    i<-1
    sigmivector<-vector()
    
    for (i in 1:r)
    {
      sigmi<-sum(taxacppp[i,]*taxacppp[r+1,])/sum(taxacppp[i,])
      sigmivector<-append(sigmivector,sigmi)
      i<-i+1
    }
    
    # calculating the c-p-structure
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    rmean<-nrow(taxamean)
    taxameancp<-rbind(taxamean,cpvector)
    taxameancp<-subset(taxameancp,select=which(taxameancp[l+1,]>0))
    
    taxameancp1<-subset(taxameancp,select=which(taxameancp[l+1,]==1))
    i<-1
    cp1percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp1percent<-100*sum(as.numeric(taxameancp1[i,]))/sum(as.numeric(taxameancp[i,]))
      cp1percentvector<-append(cp1percentvector,cp1percent)
      i<-i+1
    }
    
    taxameancp2<-subset(taxameancp,select=which(taxameancp[l+1,]==2))
    i<-1
    cp2percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp2percent<-100*sum(as.numeric(taxameancp2[i,]))/sum(as.numeric(taxameancp[i,]))
      cp2percentvector<-append(cp2percentvector,cp2percent)
      i<-i+1
    }
    
    taxameancp3<-subset(taxameancp,select=which(taxameancp[l+1,]==3))
    i<-1
    cp3percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp3percent<-100*sum(as.numeric(taxameancp3[i,]))/sum(as.numeric(taxameancp[i,]))
      cp3percentvector<-append(cp3percentvector,cp3percent)
      i<-i+1
    }
    
    taxameancp4<-subset(taxameancp,select=which(taxameancp[l+1,]==4))
    i<-1
    cp4percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp4percent<-100*sum(as.numeric(taxameancp4[i,]))/sum(as.numeric(taxameancp[i,]))
      cp4percentvector<-append(cp4percentvector,cp4percent)
      i<-i+1
    }
    
    taxameancp5<-subset(taxameancp,select=which(taxameancp[l+1,]==5))
    i<-1
    cp5percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp5percent<-100*sum(as.numeric(taxameancp5[i,]))/sum(as.numeric(taxameancp[i,]))
      cp5percentvector<-append(cp5percentvector,cp5percent)
      i<-i+1
    }
    
    # subsetting taxa which have c-p-values and the calculation of MI
    
    taxami<-taxacp
    
    i<-1
    mivector<-vector()
    
    for (i in 1:r)
    {
      mi<-sum(taxami[i,]*taxami[r+1,])/sum(taxami[i,])
      mivector<-append(mivector,mi)
      i<-i+1
    }
    
    # subsetting those taxa having c-p-values 2-5 and the calculation of MI2-5
    
    r<-nrow(taxa)
    
    taxami25<-subset(taxacp,select=which(taxacp[r+1,]>1))
    
    i<-1
    mi25vector<-vector()
    
    for (i in 1:r)
    {
      mi25<-sum(taxami25[i,]*taxami25[r+1,])/sum(taxami25[i,])
      mi25vector<-append(mi25vector,mi25)
      i<-i+1
    }
    
    # extracting p-p values from the database
    
    c<-ncol(taxa)
    i<-1
    ppvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        taxon<-names(taxa)[i]
        pp<-dat[dat["taxon"]==taxon][2]
      } else {
        pp<-0
      }
      ppvector<-append(ppvector,pp)
    }
    
    ppvector<-as.numeric(ppvector)
    taxapp<-rbind(taxa,ppvector)
    taxapp<-subset(taxapp,select=which(taxapp[r+1,]>0))
    
    # calculating the p-p-structure
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    rmean<-nrow(taxamean)
    taxameanpp<-rbind(taxamean,ppvector)
    taxameanpp<-subset(taxameanpp,select=which(taxameanpp[l+1,]>0))
    
    taxameanpp2<-subset(taxameanpp,select=which(taxameanpp[l+1,]==2))
    i<-1
    pp2percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp2percent<-100*sum(as.numeric(taxameanpp2[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp2percentvector<-append(pp2percentvector,pp2percent)
      i<-i+1
    }
    
    taxameanpp3<-subset(taxameanpp,select=which(taxameanpp[l+1,]==3))
    i<-1
    pp3percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp3percent<-100*sum(as.numeric(taxameanpp3[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp3percentvector<-append(pp3percentvector,pp3percent)
      i<-i+1
    }
    
    taxameanpp4<-subset(taxameanpp,select=which(taxameanpp[l+1,]==4))
    i<-1
    pp4percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp4percent<-100*sum(as.numeric(taxameanpp4[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp4percentvector<-append(pp4percentvector,pp4percent)
      i<-i+1
    }
    
    taxameanpp5<-subset(taxameanpp,select=which(taxameanpp[l+1,]==5))
    i<-1
    pp5percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp5percent<-100*sum(as.numeric(taxameanpp5[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp5percentvector<-append(pp5percentvector,pp5percent)
      i<-i+1
    }
    
    # subsetting those taxa having p-p-values and the calculation of PPI
    
    r<-nrow(taxa)
    
    taxappi<-taxapp
    
    i<-1
    ppivector<-vector()
    
    for (i in 1:r)
    {
      ppi<-sum(taxappi[i,]*taxappi[r+1,])/sum(taxappi[i,])
      ppivector<-append(ppivector,ppi)
      i<-i+1
    }
    
    # the calculation of the feeding type composition
    
    taxafeed<-rbind(taxa,feedvector)
    taxafeed<-subset(taxafeed,select=which(taxafeed[r+1,]>0))
    
    pf<-which(taxafeed[r+1,]<2)
    ff<-which(taxafeed[r+1,]==2)
    bf<-which(taxafeed[r+1,]==3)
    pred<-which(taxafeed[r+1,]==5)
    ov<-which(taxafeed[r+1,]==8)
    l<-length(levels(treats))
    
    taxamean<-means.by(taxanorm,treats)
    taxameanfeed<-rbind(taxamean,feedvector,plantvector)
    taxameanfeed<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]>0))
    taxameanfeedfl<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]>=2))
    taxameanfeedh<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]==1))
    
    rmean<-nrow(taxamean)
    
    taxameanpf<-subset(taxameanfeed,select=pf)
    
    i<-1
    pfpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      pfpercent<-100*sum(as.numeric(taxameanpf[i,]))/sum(as.numeric(taxameanfeed[i,]))
      pfpercentvector<-append(pfpercentvector,pfpercent)
      i<-i+1
    }
    
    taxameanff<-subset(taxameanfeed,select=ff)
    
    i<-1
    ffpercentvector<-vector()
    ffpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      ffpercent<-100*sum(as.numeric(taxameanff[i,]))/sum(as.numeric(taxameanfeed[i,]))
      ffpercentfl<-100*sum(as.numeric(taxameanff[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      ffpercentvector<-append(ffpercentvector,ffpercent)
      ffpercentvectorfl<-append(ffpercentvectorfl,ffpercentfl)
      i<-i+1
    }
    
    taxameanbf<-subset(taxameanfeed,select=bf)
    
    i<-1
    bfpercentvector<-vector()
    bfpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      bfpercent<-100*sum(as.numeric(taxameanbf[i,]))/sum(as.numeric(taxameanfeed[i,]))
      bfpercentfl<-100*sum(as.numeric(taxameanbf[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      bfpercentvector<-append(bfpercentvector,bfpercent)
      bfpercentvectorfl<-append(bfpercentvectorfl,bfpercentfl)
      i<-i+1
    }
    
    taxameanpred<-subset(taxameanfeed,select=pred)
    
    i<-1
    predpercentvector<-vector()
    predpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      predpercent<-100*sum(as.numeric(taxameanpred[i,]))/sum(as.numeric(taxameanfeed[i,]))
      predpercentfl<-100*sum(as.numeric(taxameanpred[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      predpercentvector<-append(predpercentvector,predpercent)
      predpercentvectorfl<-append(predpercentvectorfl,predpercentfl)
      i<-i+1
    }
    
    taxameanov<-subset(taxameanfeed,select=ov)
    
    i<-1
    ovpercentvector<-vector()
    ovpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      ovpercent<-100*sum(as.numeric(taxameanov[i,]))/sum(as.numeric(taxameanfeed[i,]))
      ovpercentfl<-100*sum(as.numeric(taxameanov[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      ovpercentvector<-append(ovpercentvector,ovpercent)
      ovpercentvectorfl<-append(ovpercentvectorfl,ovpercentfl)
      i<-i+1
    }
    
    # the calculation of the herbivore feeding type composition
    
    a<-which(taxameanfeedh[l+2,]=="a")
    b<-which(taxameanfeedh[l+2,]=="b")
    c<-which(taxameanfeedh[l+2,]=="c")
    d<-which(taxameanfeedh[l+2,]=="d")
    e<-which(taxameanfeedh[l+2,]=="e")
    f<-which(taxameanfeedh[l+2,]=="f")
    
    taxameana<-subset(taxameanfeedh,select=a)
    
    i<-1
    apercentvector<-vector()
    
    for (i in 1:rmean)
    {
      apercent<-100*sum(as.numeric(taxameana[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      apercentvector<-append(apercentvector,apercent)
      i<-i+1
    }
    
    taxameanb<-subset(taxameanfeedh,select=b)
    
    i<-1
    bpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      bpercent<-100*sum(as.numeric(taxameanb[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      bpercentvector<-append(bpercentvector,bpercent)
      i<-i+1
    }
    
    taxameanc<-subset(taxameanfeedh,select=c)
    
    i<-1
    cpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      cpercent<-100*sum(as.numeric(taxameanc[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      cpercentvector<-append(cpercentvector,cpercent)
      i<-i+1
    }
    
    taxameand<-subset(taxameanfeedh,select=d)
    
    i<-1
    dpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      dpercent<-100*sum(as.numeric(taxameand[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      dpercentvector<-append(dpercentvector,dpercent)
      i<-i+1
    }
    
    taxameane<-subset(taxameanfeedh,select=e)
    
    i<-1
    epercentvector<-vector()
    
    for (i in 1:rmean)
    {
      epercent<-100*sum(as.numeric(taxameane[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      epercentvector<-append(epercentvector,epercent)
      i<-i+1
    }
    
    taxameanf<-subset(taxameanfeedh,select=f)
    
    i<-1
    fpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      fpercent<-100*sum(as.numeric(taxameanf[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      fpercentvector<-append(fpercentvector,fpercent)
      i<-i+1
    }
    
    # extracting food web data from the input
    
    feedvectorcopy<-feedvector
    feedvectorcopy[which(feedvector==1)]<-"PF"
    feedvectorcopy[which(feedvector==2)]<-"FF"
    feedvectorcopy[which(feedvector==3)]<-"BF"
    feedvectorcopy[which(feedvector==5)]<-"PRED"
    feedvectorcopy[which(feedvector==7)]<-"PAR"
    feedvectorcopy[which(feedvector==8)]<-"OV"
    groupvector<-paste(feedvectorcopy,cpvector,sep="")
    
    c<-ncol(taxa)
    i<-1
    fwdvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        fwd<-0
      } else {
        if (cpvector[i]==1 || cpvector[i]==4) {
          fwd<-3.2
        }
        if (cpvector[i]==2) {
          fwd<-0.8
        }
        
        if (cpvector[i]==3) {
          fwd<-1.8
        }
        if (cpvector[i]==5) {
          fwd<-5.0
        }
      }
      fwdvector<-append(fwdvector,fwd)
    }
    
    i<-1
    footvector<-vector()
    massvector<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      mass<-as.numeric(dat[dat["taxon"]==taxon][4])
      massvector<-append(massvector,mass)
      cp<-as.numeric(dat[dat["taxon"]==taxon][2])
      foot<-0.1*(mass/cp)+0.273*(mass^0.75)
      footvector<-append(footvector,foot)
      i<-i+1
    }
    
    fwdvector<-as.numeric(fwdvector)
    taxafw<-rbind(taxa,groupvector,fwdvector,footvector,massvector)
    
    # calculation of the composite footprint and total biomass
    
    r<-nrow(taxa)
    
    i<-1
    footcvector<-vector()
    totalmassvector<-vector()
    
    for (i in 1:r)
    {
      footc<-sum(as.numeric(taxafw[i,])*as.numeric(taxafw[r+3,]))
      footcvector<-append(footcvector,footc)
      totalmass<-sum(as.numeric(taxafw[i,])*as.numeric(taxafw[r+4,]))/1000
      totalmassvector<-append(totalmassvector,totalmass)
      i<-i+1
    }
    
    # calculation of the total numbers
    
    i<-1
    totalnumvector<-vector()
    
    for (i in 1:r)
    {
      totalnum<-sum(as.numeric(taxafw[i,]))
      totalnumvector<-append(totalnumvector,totalnum)
      i<-i+1
    }
    
    # calculation of the herbivore footprint
    
    taxah<-subset(taxafw,select=which(taxafw[r+2,]==0))
    
    r<-nrow(taxa)
    
    i<-1
    foothvector<-vector()
    
    for (i in 1:r)
    {
      footh<-sum(as.numeric(taxah[i,])*as.numeric(taxah[r+3,]))
      foothvector<-append(foothvector,footh)
      i<-i+1
    }
    
    # subsetting freeliving taxa and the calculation of channel index
    
    taxafw<-subset(taxafw,select=which(taxafw[r+2,]>0))
    
    pred1<-which(taxafw[r+1,]=="PRED1")
    pred2<-which(taxafw[r+1,]=="PRED2")
    pred3<-which(taxafw[r+1,]=="PRED3")
    pred4<-which(taxafw[r+1,]=="PRED4")
    pred5<-which(taxafw[r+1,]=="PRED5")
    ov1<-which(taxafw[r+1,]=="OV1")
    ov2<-which(taxafw[r+1,]=="OV2")
    ov3<-which(taxafw[r+1,]=="OV3")
    ov4<-which(taxafw[r+1,]=="OV4")
    ov5<-which(taxafw[r+1,]=="OV5")
    ff1<-which(taxafw[r+1,]=="FF1")
    ff2<-which(taxafw[r+1,]=="FF2")
    ff3<-which(taxafw[r+1,]=="FF3")
    ff4<-which(taxafw[r+1,]=="FF4")
    ff5<-which(taxafw[r+1,]=="FF5")
    bf1<-which(taxafw[r+1,]=="BF1")
    bf2<-which(taxafw[r+1,]=="BF2")
    bf3<-which(taxafw[r+1,]=="BF3")
    bf4<-which(taxafw[r+1,]=="BF4")
    bf5<-which(taxafw[r+1,]=="BF5")
    pf1<-which(taxafw[r+1,]=="PF1")
    pf2<-which(taxafw[r+1,]=="PF2")
    pf3<-which(taxafw[r+1,]=="PF3")
    pf4<-which(taxafw[r+1,]=="PF4")
    pf5<-which(taxafw[r+1,]=="PF5")
    
    taxaff2<-subset(taxafw,select=ff2)
    taxabf1<-subset(taxafw,select=bf1)
    
    i<-1
    civector<-vector()
    
    for (i in 1:r)
    {
      ci<-100*sum(as.numeric(taxaff2[r+2,])*as.numeric(taxaff2[i,]))/(sum(as.numeric(taxaff2[r+2,])*as.numeric(taxaff2[i,]))+sum(as.numeric(taxabf1[r+2,])*as.numeric(taxabf1[i,])))
      civector<-append(civector,ci)
      i<-i+1
    }
    
    # calculation of the fungivore footprint
    
    taxaf<-subset(taxafw,select=c(ff1,ff2,ff3,ff4,ff5))
    
    r<-nrow(taxa)
    
    i<-1
    footfvector<-vector()
    
    for (i in 1:r)
    {
      footf<-sum(as.numeric(taxaf[i,])*as.numeric(taxaf[r+3,]))
      footfvector<-append(footfvector,footf)
      i<-i+1
    }
    
    # calculation of the bacterivore footprint
    
    taxab<-subset(taxafw,select=c(bf1,bf2,bf3,bf4,bf5))
    
    r<-nrow(taxa)
    
    i<-1
    footbvector<-vector()
    
    for (i in 1:r)
    {
      footb<-sum(as.numeric(taxab[i,])*as.numeric(taxab[r+3,]))
      footbvector<-append(footbvector,footb)
      i<-i+1
    }
    
    # calculation of the predator footprint
    
    taxapred<-subset(taxafw,select=c(pred1,pred2,pred3,pred4,pred5))
    
    r<-nrow(taxa)
    
    i<-1
    footpredvector<-vector()
    
    for (i in 1:r)
    {
      footpred<-sum(as.numeric(taxapred[i,])*as.numeric(taxapred[r+3,]))
      footpredvector<-append(footpredvector,footpred)
      i<-i+1
    }
    
    # calculation of the omnivore footprint
    
    taxaov<-subset(taxafw,select=c(ov1,ov2,ov3,ov4,ov5))
    
    r<-nrow(taxa)
    
    i<-1
    footovvector<-vector()
    
    for (i in 1:r)
    {
      footov<-sum(as.numeric(taxaov[i,])*as.numeric(taxaov[r+3,]))
      footovvector<-append(footovvector,footov)
      i<-i+1
    }
    
    # subsetting different food web components
    
    e<-c(bf1,ff2,pred1,ov1)
    taxae<-subset(taxafw,select=e)
    
    b<-c(bf2,ff2)
    taxab<-subset(taxafw,select=b)
    
    s<-c(pred2,pred3,pred4,pred5,ov2,ov3,ov4,ov5,ff3,ff4,ff5,bf3,bf4,bf5)
    taxas<-subset(taxafw,select=s)
    
    # calculation of the enrichment and structure footprints
    
    i<-1
    footevector<-vector()
    
    for (i in 1:r)
    {
      foote<-sum(as.numeric(taxae[i,])*as.numeric(taxae[r+3,]))
      footevector<-append(footevector,foote)
      i<-i+1
    }
    
    i<-1
    footsvector<-vector()
    
    for (i in 1:r)
    {
      foots<-sum(as.numeric(taxas[i,])*as.numeric(taxas[r+3,]))
      footsvector<-append(footsvector,foots)
      i<-i+1
    }
    
    k<-max(append(footsvector,footevector))/15
    
    # calculation of the enrichment, structure and basal indices
    
    i<-1
    bivector<-vector()
    sivector<-vector()
    eivector<-vector()
    
    for (i in 1:r)
    {
      bi<-100*sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,]))/(sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,]))+sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,])))
      si<-100*sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))/(sum(as.numeric(taxas[r+2,])*as.numeric(taxas[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      ei<-100*sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))/(sum(as.numeric(taxae[r+2,])*as.numeric(taxae[i,]))+sum(as.numeric(taxab[r+2,])*as.numeric(taxab[i,])))
      bivector<-append(bivector,bi)
      sivector<-append(sivector,si)
      eivector<-append(eivector,ei)
      i<-i+1
    }
    
    # ANOVA based on all the calculated indices
    
    if (is.null(attr(try(aov(sigmivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(sigmivector~treats))[[1]][5][1,],T),"condition"))==T) {
      sigmianova<-summary(aov(sigmivector~treats))[[1]][5][1,]
    } else {
      sigmianova<-NA
    }
    
    if (is.null(attr(try(aov(mivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(mivector~treats))[[1]][5][1,],T),"condition"))==T) {
      mianova<-summary(aov(mivector~treats))[[1]][5][1,]
    } else {
      mianova<-NA
    }
    
    if (is.null(attr(try(aov(mi25vector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(mi25vector~treats))[[1]][5][1,],T),"condition"))==T) {
      mi25anova<-summary(aov(mi25vector~treats))[[1]][5][1,]
    } else {
      mi25anova<-NA
    } 
    
    if (is.null(attr(try(aov(ppivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(ppivector~treats))[[1]][5][1,],T),"condition"))==T) {
      ppianova<-summary(aov(ppivector~treats))[[1]][5][1,]
    } else {
      ppianova<-NA
    }
    
    if (is.null(attr(try(aov(civector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(civector~treats))[[1]][5][1,],T),"condition"))==T) {
      cianova<-summary(aov(civector~treats))[[1]][5][1,]
    } else {
      cianova<-NA
    }
    
    if (is.null(attr(try(aov(bivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(bivector~treats))[[1]][5][1,],T),"condition"))==T) {
      bianova<-summary(aov(bivector~treats))[[1]][5][1,]
    } else {
      bianova<-NA
    }
    
    if (is.null(attr(try(aov(eivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(eivector~treats))[[1]][5][1,],T),"condition"))==T) {
      eianova<-summary(aov(eivector~treats))[[1]][5][1,]
    } else {
      eianova<-NA
    }
    
    if (is.null(attr(try(aov(sivector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(sivector~treats))[[1]][5][1,],T),"condition"))==T) {
      sianova<-summary(aov(sivector~treats))[[1]][5][1,]
    } else {
      sianova<-NA
    }
    
    if (is.null(attr(try(aov(totalmassvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(totalmassvector~treats))[[1]][5][1,],T),"condition"))==T) {
      totalmassanova<-summary(aov(totalmassvector~treats))[[1]][5][1,]  
    } else {
      totalmassanova<-NA
    }
    
    if (is.null(attr(try(aov(footcvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footcvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footcanova<-summary(aov(footcvector~treats))[[1]][5][1,]
    } else {
      footcanova<-NA
    }
    
    if (is.null(attr(try(aov(footevector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footevector~treats))[[1]][5][1,],T),"condition"))==T) {
      footeanova<-summary(aov(footevector~treats))[[1]][5][1,]
    } else {
      footeanova<-NA
    }
    
    if (is.null(attr(try(aov(footsvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footsvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footsanova<-summary(aov(footsvector~treats))[[1]][5][1,]
    } else {
      footsanova<-NA
    }
    
    if (is.null(attr(try(aov(foothvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(foothvector~treats))[[1]][5][1,],T),"condition"))==T) {
      foothanova<-summary(aov(foothvector~treats))[[1]][5][1,]
    } else {
      foothanova<-NA
    }
    
    if (is.null(attr(try(aov(footfvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footfvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footfanova<-summary(aov(footfvector~treats))[[1]][5][1,]
    } else {
      footfanova<-NA
    }
    
    if (is.null(attr(try(aov(footbvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footbvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footbanova<-summary(aov(footbvector~treats))[[1]][5][1,]
    } else {
      footbanova<-NA
    }
    
    if (is.null(attr(try(aov(footpredvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footpredvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footpredanova<-summary(aov(footpredvector~treats))[[1]][5][1,]
    } else {
      footpredanova<-NA
    }
    
    if (is.null(attr(try(aov(footovvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(footovvector~treats))[[1]][5][1,],T),"condition"))==T) {
      footovanova<-summary(aov(footovvector~treats))[[1]][5][1,]
    } else {
      footovanova<-NA
    }
    
    if (is.null(attr(try(aov(totalnumvector~treats),T),"condition"))==T && is.null(attr(try(summary(aov(totalnumvector~treats))[[1]][5][1,],T),"condition"))==T) {
      totalnumanova<-summary(aov(totalnumvector~treats))[[1]][5][1,]
    } else {
      totalnumanova<-NA
    }
    
    anovavector<-vector()
    anovavector<-append(anovavector,values=c(mianova,mi25anova,sigmianova,ppianova,cianova,bianova,eianova,sianova,totalmassanova,footcanova,footeanova,footsanova,foothanova,footfanova,footbanova,footpredanova,footovanova,totalnumanova))
    
    # preparing output
    
    final<-matrix(cbind(mivector,mi25vector,sigmivector,ppivector,civector,bivector,eivector,sivector,totalmassvector,footcvector,footevector,footsvector,foothvector,footfvector,footbvector,footpredvector,footovvector,totalnumvector),nrow=r,ncol=18)
    finalmean<-format(round(means.by(final,treats),2),nsmall=2)
    finalsd<-format(round(sd.by(final,treats),2),nsmall=2)
    ppimean<-format(round(means.by(cbind(ppivector[which(ppivector>=2)]),treats[which(ppivector>=2)]),2),nsmall=2)
    ppisd<-format(round(sd.by(cbind(ppivector[which(ppivector>=2)]),treats[which(ppivector>=2)]),2),nsmall=2)
    final<-cbind(finalmean[,1],finalsd[,1],finalmean[,2],finalsd[,2],finalmean[,3],finalsd[,3],ppimean,ppisd,finalmean[,5],finalsd[,5],finalmean[,6],finalsd[,6],finalmean[,7],finalsd[,7],finalmean[,8],finalsd[,8],finalmean[,9],finalsd[,9],finalmean[,10],finalsd[,10],finalmean[,11],finalsd[,11],finalmean[,12],finalsd[,12],finalmean[,13],finalsd[,13],finalmean[,14],finalsd[,14],finalmean[,15],finalsd[,15],finalmean[,16],finalsd[,16],finalmean[,17],finalsd[,17],finalmean[,18],finalsd[,18])
    finalanova<-format(round(anovavector,3),nsmall=3)
    
    i<-1
    for (i in 1:18) {
      if (finalanova[i]=="0.000") finalanova[i]<-"<0.001" 
      i<-i+1
    }
    
    finalanova<-c(finalanova[1],"-",finalanova[2],"-",finalanova[3],"-",finalanova[4],"-",finalanova[5],"-",finalanova[6],"-",finalanova[7],"-",finalanova[8],"-",finalanova[9],"-",finalanova[10],"-",finalanova[11],"-",finalanova[12],"-",finalanova[13],"-",finalanova[14],"-",finalanova[15],"-",finalanova[16],"-",finalanova[17],"-",finalanova[18],"-")
    final<-rbind(final,finalanova)
    final<-cbind(final,append(format(round(pfpercentvector,1),nsmall=1),"-"),append(format(round(ffpercentvector,1),nsmall=1),"-"),append(format(round(ffpercentvectorfl,1),nsmall=1),"-"),append(format(round(bfpercentvector,1),nsmall=1),"-"),append(format(round(bfpercentvectorfl,1),nsmall=1),"-"),append(format(round(predpercentvector,1),nsmall=1),"-"),append(format(round(predpercentvectorfl,1),nsmall=1),"-"),append(format(round(ovpercentvector,1),nsmall=1),"-"),append(format(round(ovpercentvectorfl,1),nsmall=1),"-"),
                 append(format(round(apercentvector,1),nsmall=1),"-"),append(format(round(bpercentvector,1),nsmall=1),"-"),append(format(round(cpercentvector,1),nsmall=1),"-"),append(format(round(dpercentvector,1),nsmall=1),"-"),append(format(round(epercentvector,1),nsmall=1),"-"),append(format(round(fpercentvector,1),nsmall=1),"-"),append(format(round(cp1percentvector,1),nsmall=1),"-"),append(format(round(cp2percentvector,1),nsmall=1),"-"),append(format(round(cp3percentvector,1),nsmall=1),"-"),append(format(round(cp4percentvector,1),nsmall=1),"-"),append(format(round(cp5percentvector,1),nsmall=1),"-"),append(format(round(pp2percentvector,1),nsmall=1),"-"),append(format(round(pp3percentvector,1),nsmall=1),"-"),append(format(round(pp4percentvector,1),nsmall=1),"-"),append(format(round(pp5percentvector,1),nsmall=1),"-"))
    final<-t(final)
    rownames.final<-c("Maturity Index [MEAN]","Maturity Index [SD]","Maturity Index 2-5 [MEAN]","Maturity Index 2-5 [SD]","Sigma Maturity Index [MEAN]","Sigma Maturity Index [SD]","Plant Parasitic Index [MEAN]","Plant Parasitic Index [SD]","Channel Index [MEAN]","Channel Index [SD]","Basal Index [MEAN]","Basal Index [SD]","Enrichment Index [MEAN]","Enrichment Index [SD]","Structure Index [MEAN]","Structure Index [SD]","Total biomass, mg [MEAN]","Total biomass, mg [SD]","Composite footprint [MEAN]","Composite footprint [SD]","Enrichment footprint [MEAN]","Enrichment footprint [SD]","Structure footprint [MEAN]","Structure footprint [SD]","Herbivore footprint [MEAN]","Herbivore footprint [SD]","Fungivore footprint [MEAN]","Fungivore footprint [SD]","Bacterivore footprint [MEAN]","Bacterivore footprint [SD]","Predator footprint [MEAN]","Predator footprint [SD]","Omnivore footprint [MEAN]","Omnivore footprint [SD]","Total number, ind [MEAN]","Total number, ind [SD]","Herbivores, % of total","Fungivores, % of total","Fungivores, % of free-living","Bacterivores, % of total","Bacterivores, % of free-living","Predators, % of total","Predators, % of free-living","Omnivores, % of total","Omnivores, % of free-living",
                      "Sedentary parasites, % of herbivores","Migratory endoparasites, % of herbivores","Semi-endoparasites, % of herbivores","Ectoparasites, % of herbivores","Epidermal/root hair feeders, % of herbivores","Algal/lichen/moss feeders, % of herbivores","CP 1, % of free-living","CP 2, % of free-living","CP 3, % of free-living","CP 4, % of free-living","CP 5, % of free-living","PP 2, % of herbivores","PP 3, % of herbivores","PP 4, % of herbivores","PP 5, % of herbivores")
    final<-cbind(rownames.final,final)
    colnames(final)<-c("Index name",as.character(levels(treats)),"ANOVA,p")
    
    final
  },include.rownames=F)
  
  
  
  ### [START]
  
  
  
  output$warnings<-renderTable(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Checking taxa names for repeats and spelling mistakes',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    duplindicator<-0
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    
    taxa<-tax[,-1]
    treats<-as.factor(tax[,1])
    
    c<-ncol(taxa)
    i<-1
    found<-vector()
    guess<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      i<-i+1
      cp<-dat[dat["taxon"]==taxon][2]
      if (is.na(cp)==TRUE) {
        found.value<-taxon
        found<-append(found,found.value)
        guess.value<-as.character(dat[,1][agrep(taxon,dat[,1])])[1]
        guess<-append(guess,guess.value)
      }
    }
    
    
    
    if (length(found)==0) {
      
      i<-1
      cyst<-0
      
      for (i in 1:c)
      {
        taxon<-names(taxa)[i]
        if (dat[dat["taxon"]==taxon][5]=="a") {
          cyst<-1
        }
      }
      
      if (duplindicator==0 && length(grep("dauer",names(taxa)))==0) {
        nowarnings<-matrix(ncol=1,nrow=0)  
      } 
      if (duplindicator==0 && length(grep("dauer",names(taxa)))>0) {
        nowarnings<-matrix("Your data set includes dauer larvae. You can optionally exclude them from calculations ticking a box in the sidebar panel.",ncol=1,nrow=1)  
      }
      if (duplindicator==1 && length(grep("dauer",names(taxa)))==0) {
        nowarnings<-matrix("There are some columns with identical taxa names. They have been summed.",ncol=1,nrow=1)  
      }
      if (duplindicator==1 && length(grep("dauer",names(taxa)))>0) {
        nowarnings<-matrix(c("There are some columns with identical taxa names. They have been summed.",
                             "Your data set includes dauer larvae. You can optionally exclude them from calculations ticking a box in the sidebar panel."),ncol=1,nrow=2)  
      }
      colnames(nowarnings)<-"All taxa names have been read correctly."
      
      if (cyst==1) {
        nowarnings<-rbind(nowarnings,
                          "Your data set contains cyst or root-knot nematodes. Mass estimates of all nematode species are based on measurements of adults but the mentioned groups have swollen females which have particularly high mass in comparison to other nematodes. Please consider this when comparing metabolic footprint values across samples.",
                          "You may proceed to the other tabs.")
      } else {
        nowarnings<-rbind(nowarnings,
                          "You may proceed to the other tabs.")
      }
      
      nowarnings
      
    } else {
      warnings<-cbind(found,guess)
      warnings<-rbind(warnings,c("PLEASE CORRECT YOUR INPUT FILE","AND UPLOAD IT AGAIN!"))
      colnames(warnings)<-c("No such name(s) in the database:","Suggestion for correction:")
      warnings
    }
    
  })
  
  
  
  ### [Side panel]:"HOW TO FILL"
  
  
  
  output$howtofill<-renderText(function(){
    return('<img src="howtofill.png" alt="How to fill">') 
  })
  
  
  
  ### [Feeding types & c-p/p-p]:Table with the database info
  
  
  
  output$list<-renderTable(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Extracting database information',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    
    # extracting feeding types from the database and sorting
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      if (dat[dat["taxon"]==taxon][5]=="") {
        feed<-dat[dat["taxon"]==taxon][3]
      } else {
        feed<-paste(dat[dat["taxon"]==taxon][3],dat[dat["taxon"]==taxon][5],sep="")
      }
      feedvector<-append(feedvector,feed)
      i<-i+1
    }
    
    sort<-order(substr(feedvector,1,1),names(taxa))
    names(taxa)<-names(taxa)[sort]
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in sort)
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # extracting p-p values from the database
    
    i<-1
    ppvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        taxon<-names(taxa)[i]
        pp<-dat[dat["taxon"]==taxon][2]
      } else {
        pp<-0
      }
      ppvector<-append(ppvector,pp)
    }
    
    ppvector<-as.numeric(ppvector)
    
    #extracting mass data from the database
    
    i<-1
    massvector<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      mass<-as.numeric(dat[dat["taxon"]==taxon][4])
      massvector<-append(massvector,mass)
      i<-i+1
    }
    
    # output
    
    cpvector<-round(cpvector,0)
    ppvector<-round(ppvector,0)
    
    feedvector<-paste(feedvector,plantvector,sep="")
    
    i<-1
    for (i in 1:c) {
      if (feedvector[i]==1) feedvector[i]<-"Herbivores"
      if (feedvector[i]=="1a") feedvector[i]<-"Herbivores - sedentary parasites"
      if (feedvector[i]=="1b") feedvector[i]<-"Herbivores - migratory endoparasites"
      if (feedvector[i]=="1c") feedvector[i]<-"Herbivores - semi-endoparasites"
      if (feedvector[i]=="1d") feedvector[i]<-"Herbivores - ectoparasites"
      if (feedvector[i]=="1e") feedvector[i]<-"Herbivores - epidermal/root hair feeders"
      if (feedvector[i]=="1f") feedvector[i]<-"Herbivores - algal/lichen/moss feeders"
      if (feedvector[i]==2) feedvector[i]<-"Fungivores"
      if (feedvector[i]==3) feedvector[i]<-"Bacterivores"
      if (feedvector[i]==5) feedvector[i]<-"Predators"
      if (feedvector[i]==7) feedvector[i]<-"Animal parasites (dispersal stages)"
      if (feedvector[i]==8) feedvector[i]<-"Omnivores"
    }
    
    massvector<-format(round(massvector,3),nsmall=3)
    final<-matrix(cbind(cpvector,ppvector,feedvector,massvector),nrow=c,ncol=4)
    rownames(final)<-colnames(taxa)
    colnames(final)<-c("C-p class","P-p class","Feeding type","Mass, ug")
    
    final
  })
  
  
  
  ### [Feeding types & c-p/p-p]:Feeding types, all nematodes
  
  
  
  output$feedingall<-renderGvis(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Depicting feeding type composition 1/3',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # creating a matrix with mean values per treatment
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    taxameanfeed<-rbind(taxamean,feedvector)
    taxameanfeed<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]>0))
    
    pf<-which(taxameanfeed[l+1,]<2)
    ff<-which(taxameanfeed[l+1,]==2)
    bf<-which(taxameanfeed[l+1,]==3)
    pred<-which(taxameanfeed[l+1,]==5)
    ov<-which(taxameanfeed[l+1,]==8)
    
    # plotting the feeding type composition
    
    rmean<-nrow(taxamean)
    
    taxameanpf<-subset(taxameanfeed,select=pf)
    
    i<-1
    pfpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      pfpercent<-100*sum(as.numeric(taxameanpf[i,]))/sum(as.numeric(taxameanfeed[i,]))
      pfpercentvector<-append(pfpercentvector,pfpercent)
      i<-i+1
    }
    
    taxameanff<-subset(taxameanfeed,select=ff)
    
    i<-1
    ffpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      ffpercent<-100*sum(as.numeric(taxameanff[i,]))/sum(as.numeric(taxameanfeed[i,]))
      ffpercentvector<-append(ffpercentvector,ffpercent)
      i<-i+1
    }
    
    taxameanbf<-subset(taxameanfeed,select=bf)
    
    i<-1
    bfpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      bfpercent<-100*sum(as.numeric(taxameanbf[i,]))/sum(as.numeric(taxameanfeed[i,]))
      bfpercentvector<-append(bfpercentvector,bfpercent)
      i<-i+1
    }
    
    taxameanpred<-subset(taxameanfeed,select=pred)
    
    i<-1
    predpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      predpercent<-100*sum(as.numeric(taxameanpred[i,]))/sum(as.numeric(taxameanfeed[i,]))
      predpercentvector<-append(predpercentvector,predpercent)
      i<-i+1
    }
    
    taxameanov<-subset(taxameanfeed,select=ov)
    
    i<-1
    ovpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      ovpercent<-100*sum(as.numeric(taxameanov[i,]))/sum(as.numeric(taxameanfeed[i,]))
      ovpercentvector<-append(ovpercentvector,ovpercent)
      i<-i+1
    }
    
    final<-data.frame(as.character(levels(treats)),round(pfpercentvector,2),round(ffpercentvector,2),round(bfpercentvector,2),round(predpercentvector,2),round(ovpercentvector,2))
    colnames(final)<-c("treats","Herbivores","Fungivores","Bacterivores","Predators","Omnivores")
    
    gvisColumnChart(final,xvar="treats",yvar=c("Herbivores","Fungivores","Bacterivores","Predators","Omnivores"), options=list(isStacked=TRUE,vAxis="{title:'Fraction, % of total number'}",colors="['green','orange','red','grey','blue']",title="Feeding type composition of nematode assemblage",chartArea="{width:'50%'}"))
    
  })
  
  
  
  ### [Feeding types & c-p/p-p]:Feeding types, free-living nematodes
  
  
  
  output$feedingfl<-renderGvis(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Depicting feeding type composition 2/3',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # creating a matrix with mean values per treatment
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    taxameanfeed<-rbind(taxamean,feedvector)
    taxameanfeedfl<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]>=2))
    
    ff<-which(taxameanfeedfl[l+1,]==2)
    bf<-which(taxameanfeedfl[l+1,]==3)
    pred<-which(taxameanfeedfl[l+1,]==5)
    ov<-which(taxameanfeedfl[l+1,]==8)
    
    # plotting the feeding type composition
    
    rmean<-nrow(taxamean)
    
    taxameanff<-subset(taxameanfeedfl,select=ff)
    
    i<-1
    ffpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      ffpercentfl<-100*sum(as.numeric(taxameanff[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      ffpercentvectorfl<-append(ffpercentvectorfl,ffpercentfl)
      i<-i+1
    }
    
    taxameanbf<-subset(taxameanfeedfl,select=bf)
    
    i<-1
    bfpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      bfpercentfl<-100*sum(as.numeric(taxameanbf[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      bfpercentvectorfl<-append(bfpercentvectorfl,bfpercentfl)
      i<-i+1
    }
    
    taxameanpred<-subset(taxameanfeedfl,select=pred)
    
    i<-1
    predpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      predpercentfl<-100*sum(as.numeric(taxameanpred[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      predpercentvectorfl<-append(predpercentvectorfl,predpercentfl)
      i<-i+1
    }
    
    taxameanov<-subset(taxameanfeedfl,select=ov)
    
    i<-1
    ovpercentvectorfl<-vector()
    
    for (i in 1:rmean)
    {
      ovpercentfl<-100*sum(as.numeric(taxameanov[i,]))/sum(as.numeric(taxameanfeedfl[i,]))
      ovpercentvectorfl<-append(ovpercentvectorfl,ovpercentfl)
      i<-i+1
    }
    
    final<-data.frame(as.character(levels(treats)),round(ffpercentvectorfl,2),round(bfpercentvectorfl,2),round(predpercentvectorfl,2),round(ovpercentvectorfl,2))
    colnames(final)<-c("treats","Fungivores","Bacterivores","Predators","Omnivores")
    
    gvisColumnChart(final,xvar="treats",yvar=c("Fungivores","Bacterivores","Predators","Omnivores"), options=list(isStacked=TRUE,vAxis="{title:'Fraction, % of free-living nematodes'}",colors="['orange','red','grey','blue']",title="Feeding type composition of the free-living nematode assemblage",chartArea="{width:'50%'}"))
    
  })
  
  
  
  ### [Feeding types & c-p/p-p]:Feeding types, herbivores
  
  
  
  output$feedingh<-renderGvis(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Depicting feeding type composition 3/3',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # creating a matrix with mean values per treatment
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    taxameanfeed<-rbind(taxamean,feedvector,plantvector)
    taxameanfeedh<-subset(taxameanfeed,select=which(taxameanfeed[l+1,]==1))
    
    a<-which(taxameanfeedh[l+2,]=="a")
    b<-which(taxameanfeedh[l+2,]=="b")
    c<-which(taxameanfeedh[l+2,]=="c")
    d<-which(taxameanfeedh[l+2,]=="d")
    e<-which(taxameanfeedh[l+2,]=="e")
    f<-which(taxameanfeedh[l+2,]=="f")
    
    # plotting the feeding type composition
    
    rmean<-nrow(taxamean)
    
    taxameana<-subset(taxameanfeedh,select=a)
    
    i<-1
    apercentvector<-vector()
    
    for (i in 1:rmean)
    {
      apercent<-100*sum(as.numeric(taxameana[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      apercentvector<-append(apercentvector,apercent)
      i<-i+1
    }
    
    taxameanb<-subset(taxameanfeedh,select=b)
    
    i<-1
    bpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      bpercent<-100*sum(as.numeric(taxameanb[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      bpercentvector<-append(bpercentvector,bpercent)
      i<-i+1
    }
    
    taxameanc<-subset(taxameanfeedh,select=c)
    
    i<-1
    cpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      cpercent<-100*sum(as.numeric(taxameanc[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      cpercentvector<-append(cpercentvector,cpercent)
      i<-i+1
    }
    
    taxameand<-subset(taxameanfeedh,select=d)
    
    i<-1
    dpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      dpercent<-100*sum(as.numeric(taxameand[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      dpercentvector<-append(dpercentvector,dpercent)
      i<-i+1
    }
    
    taxameane<-subset(taxameanfeedh,select=e)
    
    i<-1
    epercentvector<-vector()
    
    for (i in 1:rmean)
    {
      epercent<-100*sum(as.numeric(taxameane[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      epercentvector<-append(epercentvector,epercent)
      i<-i+1
    }
    
    taxameanf<-subset(taxameanfeedh,select=f)
    
    i<-1
    fpercentvector<-vector()
    
    for (i in 1:rmean)
    {
      fpercent<-100*sum(as.numeric(taxameanf[i,]))/sum(as.numeric(taxameanfeedh[i,]))
      fpercentvector<-append(fpercentvector,fpercent)
      i<-i+1
    }
    
    final<-data.frame(as.character(levels(treats)),round(apercentvector,2),round(bpercentvector,2),round(cpercentvector,2),round(dpercentvector,2),round(epercentvector,2),round(fpercentvector,2))
    colnames(final)<-c("treats","Sedentary parasites","Migratory endoparasites","Semi- endoparasites","Ectoparasites","Epidermal/root hair feeders","Algal/lichen/moss feeders")
    
    gvisColumnChart(final,xvar="treats",yvar=c("Sedentary parasites","Migratory endoparasites","Semi- endoparasites","Ectoparasites","Epidermal/root hair feeders","Algal/lichen/moss feeders"), options=list(isStacked=TRUE,vAxis="{title:'Fraction, % of herbivore nematodes'}",colors="['blue','orange','red','grey','lime','brown']",title="Feeding type composition of the herbivore nematode assemblage",chartArea="{width:'50%'}"))
    
  })
  
  
  
  ### [Feeding types & c-p/p-p]:c-p-structure
  
  
  
  output$cpcomp<-renderGvis(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Depicting c-p composition',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting c-p values from the database and subsetting those having c-p
    
    c<-ncol(taxa)
    i<-1
    cpvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]>1) {
        taxon<-names(taxa)[i]
        cp<-dat[dat["taxon"]==taxon][2]
      } else {
        cp<-0
      }
      cpvector<-append(cpvector,cp)
    }
    
    cpvector<-as.numeric(cpvector)
    
    # representing the cp-structure
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    rmean<-nrow(taxamean)
    taxameancp<-rbind(taxamean,cpvector)
    taxameancp<-subset(taxameancp,select=which(taxameancp[l+1,]>0))
    
    taxameancp1<-subset(taxameancp,select=which(taxameancp[l+1,]==1))
    i<-1
    cp1percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp1percent<-100*sum(as.numeric(taxameancp1[i,]))/sum(as.numeric(taxameancp[i,]))
      cp1percentvector<-append(cp1percentvector,cp1percent)
      i<-i+1
    }
    
    taxameancp2<-subset(taxameancp,select=which(taxameancp[l+1,]==2))
    i<-1
    cp2percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp2percent<-100*sum(as.numeric(taxameancp2[i,]))/sum(as.numeric(taxameancp[i,]))
      cp2percentvector<-append(cp2percentvector,cp2percent)
      i<-i+1
    }
    
    taxameancp3<-subset(taxameancp,select=which(taxameancp[l+1,]==3))
    i<-1
    cp3percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp3percent<-100*sum(as.numeric(taxameancp3[i,]))/sum(as.numeric(taxameancp[i,]))
      cp3percentvector<-append(cp3percentvector,cp3percent)
      i<-i+1
    }
    
    taxameancp4<-subset(taxameancp,select=which(taxameancp[l+1,]==4))
    i<-1
    cp4percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp4percent<-100*sum(as.numeric(taxameancp4[i,]))/sum(as.numeric(taxameancp[i,]))
      cp4percentvector<-append(cp4percentvector,cp4percent)
      i<-i+1
    }
    
    taxameancp5<-subset(taxameancp,select=which(taxameancp[l+1,]==5))
    i<-1
    cp5percentvector<-vector()
    
    for (i in 1:rmean)
    {
      cp5percent<-100*sum(as.numeric(taxameancp5[i,]))/sum(as.numeric(taxameancp[i,]))
      cp5percentvector<-append(cp5percentvector,cp5percent)
      i<-i+1
    }
    
    final<-data.frame(as.character(levels(treats)),round(cp1percentvector,2),round(cp2percentvector,2),round(cp3percentvector,2),round(cp4percentvector,2),round(cp5percentvector,2))
    colnames(final)<-c("treats","c-p 1","c-p 2","c-p 3","c-p 4","c-p 5")
    
    gvisColumnChart(final,xvar="treats",yvar=c("c-p 1","c-p 2","c-p 3","c-p 4","c-p 5"), options=list(isStacked=TRUE,vAxis="{title:'Fraction, % of free-living nematodes'}",colors="['red','orange','grey','lime','blue']",title="Coloniser-persister structure of the free-living nematode assemblage"))
    
  })
  
  
  
  ### [Feeding types & c-p/p-p]:p-p-structure
  
  
  
  output$ppcomp<-renderGvis(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Depicting p-p composition',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    i<-1
    feedvector<-vector()
    plantvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      if (nchar(feed)==2) {
        plant<-substr(feed,2,2)
        feed<-substr(feed,1,1)
      } else {
        plant<-""
      }
      feedvector<-append(feedvector,feed)
      plantvector<-append(plantvector,plant)
      i<-i+1
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    plantvector<-c(plantvector,"b","","","","","")
    
    feedvector<-as.numeric(feedvector)
    
    # extracting p-p values from the database
    
    c<-ncol(taxa)
    i<-1
    ppvector<-vector()
    
    for (i in 1:c)
    {
      if (feedvector[i]==1) {
        taxon<-names(taxa)[i]
        pp<-dat[dat["taxon"]==taxon][2]
      } else {
        pp<-0
      }
      ppvector<-append(ppvector,pp)
    }
    
    ppvector<-as.numeric(ppvector)
    
    # representing the pp-structure
    
    i<-1
    taxanorm<-vector()
    r<-nrow(taxa)
    
    for (i in 1:r)
    {
      norm<-as.numeric(taxa[i,])/sum(as.numeric(taxa[i,]))
      taxanorm<-rbind(taxanorm,norm)
      i<-i+1
    }
    
    l<-length(levels(treats))
    taxamean<-means.by(taxanorm,treats)
    
    rmean<-nrow(taxamean)
    taxameanpp<-rbind(taxamean,ppvector)
    taxameanpp<-subset(taxameanpp,select=which(taxameanpp[l+1,]>0))
    
    taxameanpp2<-subset(taxameanpp,select=which(taxameanpp[l+1,]==2))
    i<-1
    pp2percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp2percent<-100*sum(as.numeric(taxameanpp2[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp2percentvector<-append(pp2percentvector,pp2percent)
      i<-i+1
    }
    
    taxameanpp3<-subset(taxameanpp,select=which(taxameanpp[l+1,]==3))
    i<-1
    pp3percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp3percent<-100*sum(as.numeric(taxameanpp3[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp3percentvector<-append(pp3percentvector,pp3percent)
      i<-i+1
    }
    
    taxameanpp4<-subset(taxameanpp,select=which(taxameanpp[l+1,]==4))
    i<-1
    pp4percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp4percent<-100*sum(as.numeric(taxameanpp4[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp4percentvector<-append(pp4percentvector,pp4percent)
      i<-i+1
    }
    
    taxameanpp5<-subset(taxameanpp,select=which(taxameanpp[l+1,]==5))
    i<-1
    pp5percentvector<-vector()
    
    for (i in 1:rmean)
    {
      pp5percent<-100*sum(as.numeric(taxameanpp5[i,]))/sum(as.numeric(taxameanpp[i,]))
      pp5percentvector<-append(pp5percentvector,pp5percent)
      i<-i+1
    }
    
    par(mar=c(5,4,2,12),xpd=TRUE)
    barplot(rbind(pp2percentvector,pp3percentvector,pp4percentvector,pp5percentvector),
            xlim=c(0,rmean+2), ylim=c(0,100), xlab="Treatment",ylab="Fraction, % of herbivores",main="Life strategy structure of the herbivore nematode assemblage",
            names.arg=levels(treats), legend.text=c("pp=2","pp=3","pp=4","pp=5"), 
            args.legend=list(x="right",inset=c(-0.2,0)),col=c("yellow","grey","green","blue"))
    
    final<-data.frame(as.character(levels(treats)),round(pp2percentvector,2),round(pp3percentvector,2),round(pp4percentvector,2),round(pp5percentvector,2))
    colnames(final)<-c("treats","p-p 2","p-p 3","p-p 4","p-p 5")
    
    gvisColumnChart(final,xvar="treats",yvar=c("p-p 2","p-p 3","p-p 4","p-p 5"), options=list(isStacked=TRUE,vAxis="{title:'Fraction, % of herbivore nematodes'}",colors="['orange','grey','lime','blue']",title="Life strategy structure of the herbivore nematode assemblage"))
    
  })
  
  
  
  ### [Food web diagnostics]:FWD square interpretation
  
  
  
  output$fwd<-renderText(function(){
    return('<img src="fwd.png" alt="How to interpret">')
  })
  
  
  
  ### Feeding types adjustment
  
  
  
  output$selectUI<-renderUI({ 
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Preparing optional adjustment of feeding types',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    
    c<-ncol(taxa)
    i<-1
    found<-vector()
    guess<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      i<-i+1
      feed<-dat[dat["taxon"]==taxon][3]
      if (is.na(feed)==TRUE) {
        found.value<-taxon
        found<-append(found,found.value)
      }
    }
    
    if (length(found)==0) {
      
      c<-ncol(taxa)
      i<-1
      feedvector<-vector()
      
      for (i in 1:c)
      {
        taxon<-names(taxa)[i]
        if (dat[dat["taxon"]==taxon][5]=="") {
          feed<-dat[dat["taxon"]==taxon][3]
        } else {
          feed<-paste(dat[dat["taxon"]==taxon][3],dat[dat["taxon"]==taxon][5],sep="")
        }
        feedvector<-append(feedvector,feed)
        i<-i+1
      }
      
      # sorting
      sort<-order(substr(feedvector,1,1),names(taxa))
      
      i<-1
      for (i in 1:c) {
        if (feedvector[i]==1) feedvector[i]<-"Herbivores"
        if (feedvector[i]=="1a") feedvector[i]<-"Herbivores - sedentary parasites"
        if (feedvector[i]=="1b") feedvector[i]<-"Herbivores - migratory endoparasites"
        if (feedvector[i]=="1c") feedvector[i]<-"Herbivores - semi-endoparasites"
        if (feedvector[i]=="1d") feedvector[i]<-"Herbivores - ectoparasites"
        if (feedvector[i]=="1e") feedvector[i]<-"Herbivores - epidermal/root hair feeders"
        if (feedvector[i]=="1f") feedvector[i]<-"Herbivores - algal/lichen/moss feeders"
        if (feedvector[i]==2) feedvector[i]<-"Fungivores"
        if (feedvector[i]==3) feedvector[i]<-"Bacterivores"
        if (feedvector[i]==5) feedvector[i]<-"Predators"
        if (feedvector[i]==7) feedvector[i]<-"Animal parasites (dispersal stages)"
        if (feedvector[i]==8) feedvector[i]<-"Omnivores"
      }
      
      x<-length(feedvector)
      outputlist<-list()
      i<-1
      
      for (i in sort) {
        output<-selectInput(paste("feed",i,sep=""), colnames(taxa)[i], 
                            c("Herbivores - sedentary parasites"="1a","Herbivores - migratory endoparasites"="1b","Herbivores - semi-endoparasites"="1c","Herbivores - ectoparasites"="1d","Herbivores - epidermal/root hair feeders"="1e","Herbivores - algal/lichen/moss feeders"="1f","Fungivores"=2,"Bacterivores"=3,"Predators"=5,"Animal parasites (dispersal stages)"=7,"Omnivores"=8),
                            selected=feedvector[i])
        outputlist<-append(outputlist,output)
      }
      
      outputlist
      
    } else {return(NULL)}
    
  })
  
  
  
  ### User's file upload menu
  
  
  
  output$upload<-renderUI({ 
    if (is.null(input$taxtax)) {
      list(
        helpText("STEP 1 of 2: Arrange an input table as following:"),
        htmlOutput("howtofill"),
        fileInput("taxtax","STEP 2 of 2: Save the sheet in the usual Excel format (.xls or .xlsx, point as decimal separator) and upload it here:")
      )
    } else {
      list(
        "Upload successful. Refresh the page to upload a new file.",
        checkboxInput("dauer", "Exclude dauer larvae from calculations", FALSE)
      )
    }
    
  })
  
  
  
  ### Tabs
  
  
  
  output$tabs<-renderUI({ 
    if (is.null(input$taxtax)) {
      tabsetPanel(
        tabPanel("START",tableOutput("warnings"))
      )
    } else {
      tabsetPanel(
        tabPanel("START",tableOutput("warnings"),htmlOutput("additem")),
        tabPanel("Summary", helpText("Need a hint on how to interprete this? In the next tabs you will find short tips as well as references to the original articles."),
                 tableOutput("summary")),
        tabPanel("Feeding types & c-p/p-p",helpText("INTERPRETATION TIP: Coloniser-persister classification is based on life cycle properties. Nematodes of c-p-1 are regarded as enrichment opportunists, they have short life cycles and are often found in disturbed environments. In contrast, nematodes of c-p-5 have long life cycles and tend to inhabit stable, mature ecosystems. C-p of herbivores are called p-p. References - Bongers, T., 1990. The maturity index: an ecological measure of environmental disturbance based on nematode species composition. Oecologia, 83, pp.14–19; Bongers, T. et al., 1995. Proposed changes of c-p classification for nematodes. Russian Journal of Nematology, 3(1), pp.61–62"),
                 tableOutput("list"),htmlOutput("feedingall"),htmlOutput("feedingfl"),htmlOutput("feedingh"),htmlOutput("cpcomp"),htmlOutput("ppcomp")),
        tabPanel("MI family indices", helpText("INTERPRETATION TIP: Both MI and MI2-5 indicate the degree of maturity of an ecosystem. MI responds to disturbances in general, and MI2-5 together with c-p triangle can be used to distinguish pollution effects from eutrophication effects. PPI responds to disturbances in a more complex manner. References - Ferris, H. & Bongers, T., 2009. Indices for analysis of nematode assemblages. In M. Wilson & T. Kakouli-Duarte, eds. Nematodes as Environmental Bioindicators. Wallingford: CABI, pp. 124–145; de Goede, R., Bongers, T. & Ettema, C.H., 1993. Graphical presentation and interpretation of nematode community structure: c–p triangles. Med. Fac. Landbouww. Univ. Gent, 58(2b), pp.743–750."),
                 helpText("N.B.: in boxplots, the middle line shows median, the box represents lower and upper hinges and the whiskers depict data range."),helpText("HINT: by adjusting the width of your browser window you may adjust the width of the graphs"), plotOutput("MI"), plotOutput("MI25"), plotOutput("sigMI"), plotOutput("PPI"), imageOutput("triangle")),
        tabPanel("Food web diagnostics", helpText("INTERPRETATION TIP: Higher values of CI correspond to the higher proportion of energy transformed through the 'slow' fungal decomposition channel. EI parallels the intensity of nutrient enrichment and SI correlates with the degree of maturity of an ecosystem. For more details see the interpretation scheme below. Reference - Ferris, H., Bongers, T. & de Goede, R.G.M., 2001. A framework for soil food web diagnostics: extension of the nematode faunal analysis concept. Applied Soil Ecology, 18, pp.13–29."),
                 plotOutput("BI"),plotOutput("CI"),htmlOutput("fwd"),imageOutput("square")),
        tabPanel("Metabolic footprints", helpText("INTERPRETATION TIP: Metabolic footprints quantify the amplitude of Carbon utilisation by different food web components. The point in the middle of a rhombus represents the intersection of EI and SI and length of vertical and horizontal axes of the rhombus corresponds to the footprints of enrichment and structure components respectively. Reference - Ferris, H., 2010. Form and function: Metabolic footprints of nematodes in the soil food web. European Journal of Soil Biology, 46(2), pp.97–104."),helpText("HINT: if you have to scroll to see the whole image, open it in a new tab or window of your web browser"), 
                 imageOutput("footprint")),
        tabPanel("Compost analysis",
                 helpText("IMPORTANT CONSIDERATIONS: if the input file has any fungivore species or genera, they are first concatenated into families and these families are then counted. Some families such as Tylenchidae or Anguinidae are, by default, considered herbivorous by NINJA. However, for the calculation of compost MI (and only for this), these families are counted as fungivorous. The same refers to the genus Litylenchus - it will be considered fungivorous for this analysis, but herbivorous for the other analyses. If a taxon is considered fungivore by NINJA for this particular type of analysis (and only for this), it is NOT POSSIBLE to adjust its feeding type."),
                 helpText("Taxa from the following families are counted as 'diplogastrids': Diplogastridae, Diplogasteroididae, Neodiplogastridae, Odontopharyngidae, Tylopharyngidae."),
                 tableOutput("compostMI")),
        tabPanel("ABOUT",htmlOutput("about"))
      )
    }
    
  })
  
  
  
  ### Add items to the database
  
  output$additem<-renderUI({
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    duplindicator<-0
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    
    taxa<-tax[,-1]
    treats<-as.factor(tax[,1])
    
    c<-ncol(taxa)
    i<-1
    found<-vector()
    guess<-vector()
    
    for (i in 1:c)
    {
      taxon<-names(taxa)[i]
      i<-i+1
      cp<-dat[dat["taxon"]==taxon][2]
      if (is.na(cp)==TRUE) {
        found.value<-taxon
        found<-append(found,found.value)
        guess.value<-as.character(dat[,1][agrep(taxon,dat[,1])])[1]
        guess<-append(guess,guess.value)
      }
    }
    
    if (length(found)==0) {
      return(NULL)
    } else {
      progress <- Progress$new(session, min=1, max=10)
      on.exit(progress$close())
      
      progress$set(message = 'Preparing a form to complement the database',
                   detail = 'Please wait...')
      
      for (i in 1:10) {
        progress$set(value = i)
        Sys.sleep(0.05)
      }
      
      megainput<-function(inputId, label, data) {
        addResourcePath(
          prefix='tableinput', 
          directoryPath=system.file('tableinput', 
                                    package='shinyIncubator'))
        
        tagList(
          singleton(
            tags$head(
              tags$link(rel = 'stylesheet',
                        type = 'text/css',
                        href = 'tableinput/tableinput.css'),
              tags$script(src = 'tableinput/tableinput.js')
            )
          ),
          
          tags$div(
            class = 'control-group tableinput-container',
            tags$label(
              class = "control-label",
              label,
              tags$div(
                class = 'tableinput-buttons',
                tags$button(
                  type = 'button', class = 'btn btn-mini tableinput-settings hide',
                  tags$i(class = 'icon-cog')
                ),
                HTML('<a href="#" class="tableinput-plusrow"><i class="icon-plus-sign"></i></a>'),
                HTML('<a href="#" class="tableinput-minusrow"><i class="icon-minus-sign"></i></a>')
              )
            ),
            tags$table(
              id = inputId,
              class = 'tableinput data table table-bordered table-condensed',
              tags$colgroup(
                lapply(names(data), function(name) {
                  tags$col('data-name' = name,
                           'data-field' = name,
                           'data-type' = 'character')
                })
              ),
              tags$thead(
                class = 'hide',
                tags$tr(
                  lapply(names(data), function(name) {
                    tags$th(name)
                  })
                )
              ),
              tags$tbody(
                lapply(1:nrow(data), function(i) {
                  tags$tr(
                    lapply(names(data), function(name) {
                      tags$td(
                        div(tabindex=0, as.character(data[i,name]))
                      )
                    })
                  )
                })
              )
            ),
            tags$div(
              class = 'tableinput-editor modal hide fade',
              tags$div(
                class = 'modal-header',
                HTML('<button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>'),
                tags$h3(label)
              ),
              tags$div(
                class = 'modal-body',
                
                HTML('
                     <form class="form-horizontal">
                     <div class="control-group">
                     <label class="control-label">Rows</label>
                     <div class="controls">
                     <input type="number" class="tableinput-rowcount">
                     </div>
                     </div>
                     <div class="control-group">
                     <label class="control-label">Columns</label>
                     <div class="controls">
                     <input type="number" class="tableinput-colcount">
                     </div>
                     </div>
                     </form>'
                )
                ),
              tags$div(
                class = 'modal-footer',
                tags$a(href = '#', class = 'btn btn-primary tableinput-edit', 'OK'),
                tags$a(href = '#',
                       class = 'btn',
                       'data-dismiss' = 'modal',
                       'aria-hidden' = 'true',
                       'Cancel')
              )
                )
            )
          )
      }
      
      additem_example<-data.frame(rbind(c("Taxon name","c-p/p-p","Feeding type (as in Yeates et al. 1993)","Mass, ug","Herbivore category (as in Yeates et al. 1993)"),
                                        c("Bacterivore example","2","3","10.457",""),
                                        c("Herbivore example","3","1","5.1","b")))
      
      list(
        htmlOutput("additemexpl"),
        megainput('additem','Press (+) to add an item and type your data or (-) to remove the last item',additem_example)
      )
    }
  })
  
  
  
  output$additemexpl<-renderText({
    return(c('','<br><br>Or does the <b>database lack some taxa or contain a mistake?</b> For your current analysis, you can <b>add taxa</b> that are not contained in the NINJA database using the form below. Please note that these changes will be <b>lost</b> after your reload the page. It is therefore advised to <a href="mailto:sarasm@inia.es">report all issues</a>. We appreciate your help!<br><br>'))
  })  
  
  
  
  output$about<-renderText({
    return(c('','<center>Did you like NINJA? Do you have any suggestions for improvement?','<a href="mailto:sarasm@inia.es">Send us your feedback!</a></br></br>','If you used NINJA in your research, please refer to it by citing:</center><b>Sieriebriennikov, B., Ferris, H., and de Goede, R.G.M. (2014) <i>"NINJA: An automated calculation system for nematode-based biological monitoring".</i> European Journal of Soil Biology, 61: 90-93. <a href="http://dx.doi.org/10.1016/j.ejsobi.2014.02.004">DOI: 10.1016/j.ejsobi.2014.02.004</a></b></br></br>','Please note that the database is regularly updated (last update: Aug 29th, 2019). New taxa are added with each update, and old entries are modified if new data are published. This may affect the exact values of calculated indices!'))
  })
  
  
  
  ### Hanne's index
  
  
  
  output$compostMI<-renderTable(function(){
    if (is.null(input$taxtax)) {
      return(NULL)
    }
    
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating compost MI',
                 detail = 'Please wait...')
    
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    library(gdata)
    tax<-read.xls(input$taxtax$datapath,check.names=F)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.character)
    tax[,-(1:2)]<-sapply(tax[,-(1:2)],as.numeric)
    names(tax)[1]<-"a"
    names(tax)[2]<-"b"
    
    dup<-which(duplicated(names(tax))==T)
    
    dupl<-length(dup)
    
    namesl<-length(names(tax))
    
    if (dupl>0) {
      duplindicator<-1
      megadupset<-vector()
      for (i in 1:dupl) {
        dupset<-which(names(tax)[1:namesl]==names(tax)[dup[i]])
        megadupset<-append(megadupset,dupset)
        dupsum<-rowSums(tax[dupset])
        dupname<-names(tax)[dupset][1]
        dupsum<-cbind(dupsum)
        colnames(dupsum)<-dupname
        tax<-cbind(tax,dupsum)
      }
      
      tax<-tax[-megadupset]
      
      del<-grep(".",names(tax),fixed=T)
      
      if (length(del)>0) {
        tax<-tax[-del]
      }
    }
    
    compost.treats.ind<-as.factor(tax[,1])
    
    tax<-tax[,-1]
    dat<-read.xls("database.xls")
    if (is.null(input$additem)==F){
      datadditem<-data.frame(input$additem[-1,])
      names(datadditem)<-names(dat)
      dat<-rbind(dat,datadditem)
      dat<-dat[order(as.character(dat[,1])),]
      dat[,2]<-as.integer(dat[,2])
      dat[,3]<-as.integer(dat[,3])
      dat[,4]<-as.numeric(dat[,4])
    }
    
    taxa<-tax[,-1]
    r<-nrow(taxa)
    fakeplant<-rep(0,r)
    fakeplant<-cbind(fakeplant)
    colnames(fakeplant)<-"zzzfakeplant"
    fakefl1<-rep(0,r)
    fakefl1<-cbind(fakefl1)
    colnames(fakefl1)<-"zzzfakefl1"
    fakefl2<-rep(0,r)
    fakefl2<-cbind(fakefl2)
    colnames(fakefl2)<-"zzzfakefl2"
    fakefl3<-rep(0,r)
    fakefl3<-cbind(fakefl3)
    colnames(fakefl3)<-"zzzfakefl3"
    fakefl4<-rep(0,r)
    fakefl4<-cbind(fakefl4)
    colnames(fakefl4)<-"zzzfakefl4"
    fakefl5<-rep(0,r)
    fakefl5<-cbind(fakefl5)
    colnames(fakefl5)<-"zzzfakefl5"
    taxa<-cbind(taxa,fakeplant,fakefl1,fakefl2,fakefl3,fakefl4,fakefl5)
    treats<-as.factor(tax[,1])
    
    if (input$dauer==T && length(grep("dauer",names(taxa)))>0) {
      taxa[grep("dauer",names(taxa))]<-rep(0,r)
    } 
    if (input$dauer==T && length(grep("dauer larvae",names(taxa)))>0) {
      taxa[grep("dauer larvae",names(taxa))]<-rep(0,r)
    }
    if (input$dauer==T && length(grep("dauerlarvae",names(taxa)))>0) {
      taxa[grep("dauerlarvae",names(taxa))]<-rep(0,r)
    }
    
    # specifying some functions
    
    means.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,mean,na.rm=TRUE))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    sd.by<-function(data,INDEX){
      b<-by(data,INDEX,function(d)apply(d,2,sd))
      return(structure(
        t(matrix(unlist(b),nrow=length(b[[1]]))),
        dimnames=list(names(b),col.names=names(b[[1]]))
      ))
    }
    
    # calculation of the total numbers
    
    i<-1
    totalnumvector<-vector()
    
    for (i in 1:r)
    {
      totalnum<-sum(as.numeric(taxa[i,]))
      totalnumvector<-append(totalnumvector,totalnum)
      i<-i+1
    }
    
    # extracting feeding type data from the input
    
    c<-ncol(taxa)
    feedvector<-vector()
    
    for (i in 1:(c-6))
    {
      feed<-eval(parse(text=paste("input$feed",i,sep="")))
      feedvector<-append(feedvector,feed)
    }
    
    feedvector<-c(feedvector,1,3,2,3,5,8)
    feedvector<-as.numeric(feedvector)
    
    # preparing compost data set
    
    taxafeed<-rbind(taxa,feedvector)
    c.feed<-ncol(taxafeed)
    taxa.compost<-taxafeed[,1:(c.feed-6)]
    
    c.compost<-ncol(taxa.compost)
    ff.compost.vector<-vector()
    dipl.vector<-vector()
    
    for (i in 1:c.compost)
    {
      taxon<-names(taxa)[i]
      
      ff.compost<-dat[dat["taxon"]==taxon][7]
      ff.compost.vector<-c(ff.compost.vector,ff.compost)
      
      dipl<-dat[dat["taxon"]==taxon][8]
      dipl.vector<-c(dipl.vector,dipl)
    }
    
    taxa.compost<-rbind(taxa.compost,ff.compost.vector,dipl.vector)
    
    # counting fungivore families and identifying the presence of diplogasterids
    
    r.compost<-nrow(taxa.compost)
    n.compost.ff.fam.vector<-vector()
    dipl.presence.vector<-vector()
    
    for (i in 1:(r.compost-3))
    {
      ff.compost.i<-as.numeric(taxa.compost[i,which(ff.compost.vector!="")])
      ff.compost.i<-ff.compost.i[which(ff.compost.i>0)]
      n.compost.ff.fam<-length(ff.compost.i)
      n.compost.ff.fam.vector<-c(n.compost.ff.fam.vector,n.compost.ff.fam)
      
      dipl.i<-as.numeric(taxa.compost[i,which(dipl.vector!="")])
      dipl.i<-dipl.i[which(dipl.i>0)]
      if (length(dipl.i)>0) {
        dipl.presence<-1
      } else {
        dipl.presence<-0
      }
      dipl.presence.vector<-c(dipl.presence.vector,dipl.presence)
    }
    
    # calculating f/(f+b) ratio for compost
    
    bf.compost<-which(taxa.compost[r+1,]==3)
    ff.compost<-which(taxa.compost[r+2,]!="")
    
    r.compost<-nrow(taxa.compost)
    ff.ratio.compost.vector<-vector()
    
    for (i in 1:(r.compost-3))
    {
      n.ff.compost.i<-sum(as.numeric(taxa.compost[i,ff.compost]))
      n.bf.compost.i<-sum(as.numeric(taxa.compost[i,bf.compost]))
      if (n.bf.compost.i==0) {
        ff.ratio.compost.i<-0
      } else {
        ff.ratio.compost.i<-n.ff.compost.i/(n.ff.compost.i+n.bf.compost.i)
      }
      ff.ratio.compost.vector<-c(ff.ratio.compost.vector,ff.ratio.compost.i)
    }
    
    # calculating compost MI
    
    r.compost<-nrow(taxa.compost)
    compost.MI.vector<-vector()
    
    for (i in 1:(r.compost-3))
    {
      if (totalnumvector[i]<0.5) {
        A<-0
      }
      if (totalnumvector[i]<2000 && totalnumvector[i]>=0.5) {
        A<-1/(1+10^(-0.00381697*(totalnumvector[i]-475)))
      }
      if (totalnumvector[i]>=2000) {
        A<-1
      }
      
      if (ff.ratio.compost.vector[i]<0.05) {
        B<-0
      }
      if (ff.ratio.compost.vector[i]<0.9 && ff.ratio.compost.vector[i]>=0.05) {
        B<-1/(1+10^(-9.54243*(ff.ratio.compost.vector[i]-0.25)))
      }
      if (ff.ratio.compost.vector[i]>=0.9) {
        B<-1
      }
      
      if (n.compost.ff.fam.vector[i]==0) {
        C<-0
      }
      if (n.compost.ff.fam.vector[i]==1) {
        C<-0.3
      }
      if (n.compost.ff.fam.vector[i]==2) {
        C<-0.75
      }
      if (n.compost.ff.fam.vector[i]>=3) {
        C<-1
      }
      
      if (n.compost.ff.fam.vector[i]==0) {
        D<-0
      } else {
        D<-0.75
      }
      
      compost.MI.i<-A+B+C+D
      compost.MI.vector<-c(compost.MI.vector,compost.MI.i)
    }
    
    compost.results<-data.frame(compost.treats.ind,treats,
                                compost.MI.vector,
                                totalnumvector,ff.ratio.compost.vector,
                                n.compost.ff.fam.vector,dipl.presence.vector)
    
    colnames(compost.results)<-c("Sample","Group","Compost MI","Total number, ind","F/(F+B)","Number of fungivore families","Diplogasterids")
    a<-compost.results[,7]
    compost.results[which(a==0),7]<-"Absent"
    compost.results[which(a==1),7]<-"Present"
    
    compost.results
    
  })
  
  
  
})
