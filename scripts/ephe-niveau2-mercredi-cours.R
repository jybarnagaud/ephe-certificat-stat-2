#---------------------------------------------------------------------#
### Formation "certificat statistique pour les �cologues" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# ce script permet la r�plication des exemples du cours du mardi 
#---------------------------------------------------------------------#

rm(list=ls())
library(ggplot2)
library(cowplot)
library(ggeffects)
library(arm)
library(sjPlot)

setwd("D:/EPHE_enseignement/Statistiques/certificat/2022/niveau2/supports/donnees")

#---------------------------------------------#
### L'invasion de la processionnaire du pin ###
#---------------------------------------------#

# donn�es issues d'un suivi de l'INRA de Bordeaux et de l'ONF
chenilles=read.table("chenilles.txt",header=T,sep="\t")

# calcul de la proportion de pins attaqu�s par parcelle
chenilles$prop_attaq=chenilles$nbattaq/chenilles$nbpins
head(chenilles)

# exploration graphique
p1=ggplot(chenilles)+
		aes(x=prop_attaq)+
			geom_histogram()+
					xlab("Proportion de pins attaqu�s dans une parcelle")+
						ylab("Fr�quence")

p2=ggplot(chenilles)+
	aes(x=prop_foret,y=prop_attaq)+
		geom_point()+
			xlab("Proportion de for�t dans le paysage")+
				ylab("Proportion de pins attaqu�s dans une parcelle")

p3=ggplot(chenilles)+
	aes(x=prop_pin,y=prop_attaq)+
		geom_point()+
				xlab("Proportion de pins dans une parcelle")+
					ylab("Proportion de pins attaqu�s dans une parcelle")

p4=ggplot(chenilles)+
	aes(x=prop_foret,y=prop_pin)+
	geom_point()+
	ylab("Proportion de pins dans une parcelle")+
	xlab("Proportion de for�t dans le paysage")

cowplot::plot_grid(p1,p2,p3,p4)


# ce qu'il ne faut pas faire : mod�le lin�aire classique
mod.chenilles.faux=lm(prop_attaq~prop_foret+prop_pin,data=chenilles)
par(mfrow=c(2,2))
plot(mod.chenilles.faux)

# repr�sentation de ce mod�le faux
ef.chen.faux=effect("prop_pin",mod=mod.chenilles.faux, partial.residuals=T)
plot(ef.chen.faux, smooth.residuals=F,main="Chenille processionnaire",xlab="% de pins dans le paysage",ylab="Proportion d'arbres attaqués")

plot(ggpredict(mod.chenilles.faux,terms="prop_pin"),add.data=T)

# pr�diction fausse (probabilit� d'attaque n�gative)
predict(mod.chenilles.faux,newdata=list(prop_foret=mean(chenilles$prop_foret),prop_pin=100))

# la transformation logit
library(arm)
chenilles$logit.prop_attaq=logit(chenilles$prop_attaq)
plot(chenilles$prop_attaq,chenilles$logit.prop_attaq,xlab="Proportion d'attaques",ylab="logit(proportion d'attaques)")

# ce qu'il ne faut encore pas faire : mod�le lin�aire avec transformation logit
mod.chenilles.logit.faux=lm(logit.prop_attaq~prop_foret+prop_pin,data=chenilles)

# ce qu'il faut faire : le GLM binomial
glm.chenilles=glm(cbind(nbattaq,nbpins-nbattaq)~prop_foret+prop_pin,family=binomial,data=chenilles)

# contr�le des r�sidus
par(mfrow=c(2,2))
plot(glm.chenilles)

# interpr�tation
summary(glm.chenilles)


# proportion de d�viance expliqu�e
100*(4634.3-4235.6)/4634.3

# odds ratios (~exponentielle des param�tres)
library(questionr)
odds.ratio(glm.chenilles)

# version format�e (attention : noter l'impact de l'arrondi, l'odds ratio de la for�t est 1 et non 0.99)
tab_model(glm.chenilles)

# repr�sentation graphique sur l'�chelle naturelle, avec r�sidus partiels
p1=plot(ggemmeans(glm.chenilles,terms="prop_foret"),residuals=T)+labs(title="",x="Proportion de for�t (%)",y="Probabilit� d'infestation marginale")
p2=plot(ggemmeans(glm.chenilles,terms="prop_pin"),residuals=T)+labs(title="",x="Proportion de pins (%)",y="Probabilit� d'infestation marginale")
cowplot::plot_grid(p1,p2)

# on peut supprimer les points de r�sidus partiels afin de mieux visualiser la variation (mais attention � la range de y)
p1=plot(ggemmeans(glm.chenilles,terms="prop_foret"))+labs(title="",x="Proportion de for�t (%)",y="Probabilit� d'infestation marginale")
p2=plot(ggemmeans(glm.chenilles,terms="prop_pin"))+labs(title="",x="Proportion de pins (%)",y="Probabilit� d'infestation marginale")
cowplot::plot_grid(p1,p2)

# m�me repr�sentation avec visreg
library(visreg)
par(mfrow=c(1,2)) # sur l'�chelle du lien logit - met en valeur la relation lin�aire, mais difficile � interpr�ter
visreg(glm.chenilles,xvar="prop_foret",xlab="Proportion de for�t (%)",ylab="Logit (Probabilit� d'infestation marginale)",cex.lab=1.5)
visreg(glm.chenilles,xvar="prop_pin",xlab="Proportion de pins (%)",ylab="Logit (Probabilit� d'infestation marginale)",cex.lab=1.5)

par(mfrow=c(1,2)) # sur l'�chelle naturelle (! attention au changement d'�chelle des y)
visreg(glm.chenilles,xvar="prop_foret",scale="response",xlab="Proportion de for�t (%)",ylab="Probabilit� d'infestation marginale",cex.lab=1.5)
visreg(glm.chenilles,xvar="prop_pin",scale="response",xlab="Proportion de pins (%)",ylab="Probabilit� d'infestation marginale",cex.lab=1.5)

#------------------------------------#
### Un test de la th�orie des �les ###
#------------------------------------#

# donn�es publi�es dans Lomolino et al, Ecology, 1989
lomolino=read.table("Lomolino.txt",header=T,sep="\t")
summary(lomolino)

# exploration graphique
p1=ggplot(lomolino)+
			aes(x=DSOUR,y=NSPEC)+
				geom_point()+
					labs(title="",x="Distance � la source, km",y="Richesse sp�cifique")

p2=ggplot(lomolino)+
		aes(x=LATI,y=NSPEC)+
			geom_point()+
				labs(title="",x="Latitude",y="Richesse sp�cifique")


p3=ggplot(lomolino)+
		aes(x=SURF,y=NSPEC)+
			geom_point()+
				labs(title="",x="Surface de l'�le",y="Richesse sp�cifique")

p4=ggplot(lomolino)+
	aes(x=log(SURF),y=NSPEC)+
			geom_point()+
					labs(title="",x="log(Surface de l'�le)",y="Richesse sp�cifique")

cowplot::plot_grid(p1,p2,p3,p4)


lomolino$log.SURF=log(lomolino$SURF)

# histogramme de la variable de r�ponse
hist(lomolino$NSPEC,xlab="Nombre d'esp�ces",ylab="Fr�quence",main="")

# transfo racine carr�e
sq.NSPEC=sqrt(lomolino$NSPEC)
hist(lomolino$NSPEC,xlab="sqrt(Nombre d'esp�ces)",ylab="Fr�quence",main="")
shapiro.test(sq.NSPEC)

# transfo log
lg.NSPEC=log(lomolino$NSPEC)
hist(lg.NSPEC,xlab="log(Nombre d'esp�ces)",ylab="Fr�quence",main="")
shapiro.test(lg.NSPEC)

# mod�le avec NSPEC transform�e en racine carr�e
sq.lomo=lm(sq.NSPEC~DSOUR+LATI+log.SURF,data=lomolino)

# r�sidus de ce mod�le 
par(mfrow=c(2,2))
plot(sq.lomo)

# interpr�tation : attention, les param�tres s'entendent par rapport � la variable de r�ponse transform�e
summary(sq.lomo)

# repr�sentation graphique
plot(ggemmeans(sq.lomo,terms="DSOUR"),residuals=T)+labs(title="",x="Distance � la source (km)",y="racine(richesse sp�cifique)")
					
# la m�me avec visreg					
library(visreg)
visreg(sq.lomo,xvar="DSOUR",xlab="Distance � la source, km",ylab="sqrt(Richesse sp�cifique) marginale",cex.lab=1.5)

# pr�diction d'une richesse sp�cifique pour une distance � la source donn�e, � latitude et surface moyennes
x.pred=list(DSOUR=700,LATI=mean(lomolino$LATI),log.SURF=mean(lomolino$log.SURF))
predict(sq.lomo,newdata=x.pred) # le mod�le pr�dit une valeur impossible (racine carr�e n�gative)

# fonction de lien log
x=rnorm(1000,2,1)
plot(exp(x+1),log(x+1),xlab=expression(paste("exp(",lambda,")",sep="")),ylab=expression(paste("log(",lambda,")",sep="")),cex.lab=1.5)

# r�gression de Poisson
lomo.pois=glm(NSPEC~DSOUR+LATI+log.SURF,family="poisson",data=lomolino)

# r�sidus
par(mfrow=c(2,2))
plot(lomo.pois)


# % de d�viance expliqu�e
100*(79.8-29.1)/79.8


# loi de Poisson vs loi de Poisson surdispers�e
par(mfrow=c(1,2))
hist(rpois(1000,mean(lomolino$NSPEC)),xlab="Loi de Poisson",ylab="Fr?quence",main="")
hist(rpois(1000,mean(lomolino$NSPEC))*rbinom(1000,1,0.6),xlab="Loi de Poisson surdispers?e",ylab="Fr?quence",main="")

# surdispersion r�siduelle
res.lomo.pois=residuals(lomo.pois, type = "pearson")
N=nrow(lomolino)
k=length(coef(lomo.pois))   
sum(res.lomo.pois^2) / (N - k)

# param�tres
summary(lomo.pois)

# incidence rate ratios
library(biostat3)
eform(lomo.pois)
tab_model(lomo.pois) # attention � l'arrondi impos� par cette fonction

# displays graphiques
p1=plot(ggemmeans(lomo.pois,terms="DSOUR"),residuals=T)+labs(title="",x="Distance � la source, km",y="richesse sp�cifique marginale")
p2=plot(ggemmeans(lomo.pois,terms="log.SURF"),residuals=T)+labs(title="",x="log(surface)",y="richesse sp�cifique marginale")
cowplot::plot_grid(p1,p2)

# avec visreg, sur l'�chelle du lien et sur l'�chelle naturelle
par(mfrow=c(2,2))
visreg(lomo.pois,xvar="DSOUR",scale="linear",xlab="Distance � la source (km)",ylab="log(Richesse sp�cifique) marginale")
visreg(lomo.pois,xvar="DSOUR",scale="response",xlab="Distance � la source (km)",ylab="Richesse sp�cifique marginale")
visreg(lomo.pois,xvar="log.SURF",scale="linear",xlab="log(Surface de l'�lot, km?)",ylab="log(Richesse sp�cifique) marginale")
visreg(lomo.pois,xvar="log.SURF",scale="response",xlab="log(Surface de l'�lot, km?)",ylab="Richesse sp�cifique marginale")

# intervalles de confiance � la main (pour comprendre ce qu'il y a derri�re les codes pr�c�dents)
x.pred=list(DSOUR=seq(from=min(lomolino$DSOUR),to=max(lomolino$DSOUR),length.out=100),LATI=rep(mean(lomolino$LATI),100),log.SURF=rep(mean(lomolino$log.SURF),100))
pred.lomo.pois=predict(lomo.pois,newdata=x.pred,type="link",se.fit=T,interval="confidence")
critval <- 1.96 ## approx 95% CI

fit=pred.lomo.pois$fit
fit.lomo.pois=lomo.pois$family$linkinv(fit)

up=pred.lomo.pois$fit + (1.96 * pred.lomo.pois$se.fit)
up.lomo.pois=lomo.pois$family$linkinv(up)

low=pred.lomo.pois$fit - (1.96 * pred.lomo.pois$se.fit)
low.lomo.pois=lomo.pois$family$linkinv(low)

dev.off()
plot(x.pred$DSOUR,fit.lomo.pois,type="l",xlab="Distance � la source",ylab="Richesse sp�cifique")
lines(x.pred$DSOUR,up.lomo.pois,lty="dashed")
lines(x.pred$DSOUR,low.lomo.pois,lty="dashed")

#----------------------------------#
### En for�t avec le titipounamu ###
#----------------------------------#

rifleman=read.table("rifleman.txt",header=T,sep="\t")

# exploration des donn�es
p1=ggplot(rifleman)+
		aes(x=PROPFOR,y=total)+
			geom_point()+
				labs(title="",x="Proportion de for�t, %", y="Abondance")

p2=ggplot(rifleman)+
	aes(x=PROPNATFOR,y=total)+
		geom_point()+
			labs(title="",x="Proportion de for�t native, %", y="Abondance")

p3=ggplot(rifleman)+
	aes(x=SHA500,y=total)+
		geom_point()+
			labs(title="",x="H�t�rog�n�it� du paysage", y="Abondance")

p4=ggplot(rifleman)+
	aes(x=total)+
		geom_histogram()+
	labs(title="",x="Abondance",y="Fr�quence")

cowplot::plot_grid(p1,p2,p3,p4)

# surdispersion
mean(rifleman$total)
var(rifleman$total)
sum(rifleman$total==0)/nrow(rifleman)

# pour comparaison, loi de Poisson parfaite de m�me moyenne
rif.sim=rpois(1000,mean(rifleman$total))
hist(rif.sim,xlab="donn�es simul�es",ylab="Fr�quence",main="",cex.lab=1.5)
mean(rif.sim)
var(rif.sim)
sum(rif.sim==0)/length(rif.sim)

# glm classique
rif.pois=glm(total~PROPFOR+PROPNATFOR+SHA500,family=poisson,data=rifleman)

par(mfrow=c(2,2))
plot(rif.pois)

# zero inflation?
res.rif=residuals(rif.pois, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.pois))   
sum(res.rif^2) / (N - p)

# GLM quasi Poisson
rif.qpois=glm(total~PROPFOR+PROPNATFOR+SHA500,family=quasipoisson,data=rifleman)
par(mfrow=c(2,2))
plot(rif.qpois)
summary(rif.qpois)

res.qrif=residuals(rif.qpois, type = "pearson")
N  <- nrow(rifleman)
p  <- length(coef(rif.qpois))   
sum(res.qrif^2) / (N - p)


# GLM n�gatif binomial
library(MASS)

# une simulation pour montrer la diff�rence entre un mod�le de Poisson et une distribution n�gative binomiale
lambda.sim.pois=rpois(1000, mean(rifleman$total))
lambda.sim.nb=rnegbin(1000, mu = mean(rifleman$total), theta = 0.3)
par(mfrow=c(2,1))
hist(lambda.sim.pois,xlab="Loi de Poisson, esp�rance = 0.17",cex.lab=1.5,ylab="Fr�quence",main="")
legend("topright",legend=c(paste("moyenne =",round(mean(lambda.sim.pois),2),sep=" "),paste("variance =",round(var(lambda.sim.pois),2),sep=" ")),bty="n",cex=1.5)
hist(lambda.sim.nb,xlab="Loi N�gative Binomiale, esp�rance = 0.17, dispersion = 0.3",cex.lab=1.5,ylab="Fr�quence",main="")
legend("topright",legend=c(paste("moyenne =",round(mean(lambda.sim.nb),2),sep=" "),paste("variance =",round(var(lambda.sim.nb),2),sep=" ")),bty="n",cex=1.5)

# mod�le n�gatif binomial
rif.negbin=glm.nb(total~PROPFOR+PROPNATFOR+SHA500,data=rifleman)

# r�sidus
par(mfrow=c(2,2))
plot(rif.negbin)

# inf�rence
summary(rif.negbin)

# surdispersion
res.rif.nb=residuals(rif.negbin, type = "pearson")
sum(res.rif.nb^2) / (N - p)

# Zero Inflated Poisson
library(pscl)
rif.zip=zeroinfl(total~PROPFOR+PROPNATFOR+SHA500|region,dist="poisson",data=rifleman)

# r�sidus
plot(predict(rif.zip),residuals(rif.zip,type="pearson"),xlab="Valeurs pr�dites",ylab="Valeurs r�siduelles",cex.lab=1.5)

# surdispersion
res.rif.zip=residuals(rif.zip, type = "pearson")
p2  <- length(coef(rif.zip))   
sum(res.rif.zip^2) / (N - p2)

# inf�rence
summary(rif.zip)

# hurdle model = comme une ZIP, mais la partie poissonnienne est contrainte 
rif.hur=hurdle(total~PROPFOR+PROPNATFOR+SHA500|region,dist="poisson",data=rifleman)

# r�sidus
plot(predict(rif.hur),residuals(rif.hur,type="pearson"),xlab="Valeurs pr�dites",ylab="Valeurs r�siduelles",cex.lab=1.5)

# surdispersion
res.rif.hur=residuals(rif.hur, type = "pearson")
p3  <- length(coef(rif.hur))   
sum(res.rif.hur^2) / (N - p3)

# inf�rence
summary(rif.hur)

# repr�sentation graphique
p1=plot(ggemmeans(rif.zip,terms="PROPFOR"),show.title=F)+
xlab("Proportion de for�t")+ylab("Comptages marginaux")

p2=plot(ggemmeans(rif.zip,terms="region",type="zero_inflated"),show.title=F)+
xlab("R�gion")+ylab("Probabilit� de pr�sence marginale")

cowplot::plot_grid(p1,p2)

# comparaison des diff�rents mod�les (estimateurs et intervalles de confiance)

conf.pois=confint(rif.pois)
conf.qpois=confint(rif.qpois)
conf.negbin=confint(rif.negbin)
conf.zip=confint(rif.zip)
conf.hur=confint(rif.hur)

plot(c(1:9),c(1:9),type="n",xaxt="n",xlab="Variable",ylab="Param?tre ? IC 95%",cex.lab=1.5,ylim=c(-4,4))
axis(c(1:9),at=c(1,3,5,7,8,9),labels=c("PROPFOR","PROPNATFOR","SHA500","Foothills vs Banks","Mont vs Banks","Plain vs Banks"))
abline(h=0,lty="dashed",col="gray30")

# effet PROPFOR
arrows(x0=c(1,1.2,1.4,1.6,1.8),x1=c(1,1.2,1.4,1.6,1.8),
			 y0=c(conf.pois[2,1],conf.qpois[2,1],conf.negbin[2,1],conf.zip[2,1],conf.hur[2,1])
			 ,y1=c(conf.pois[2,2],conf.qpois[2,2],conf.negbin[2,2],conf.zip[2,2],conf.hur[2,2])
			 									,length=0,code=2,lwd=3,col=c("black","red","blue","orange","brown"))
points(c(1,1.2,1.4,1.6,1.8),c(coef(rif.pois)[2],coef(rif.qpois)[2],coef(rif.negbin)[2],coef(rif.zip)[2],coef(rif.hur)[2]),pch=21,bg=c("black","red","blue","orange","brown"),col=c("black","red","blue","orange","brown"))						

# effet PROPNATFOR
arrows(x0=c(3,3.2,3.4,3.6,3.8),x1=c(3,3.2,3.4,3.6,3.8),
			 y0=c(conf.pois[3,1],conf.qpois[3,1],conf.negbin[3,1],conf.zip[3,1],conf.hur[3,1])
			 ,y1=c(conf.pois[3,2],conf.qpois[3,2],conf.negbin[3,2],conf.zip[3,2],conf.hur[3,2])
			 ,length=0,code=2,lwd=3,col=c("black","red","blue","orange","brown"))
points(c(3,3.2,3.4,3.6,3.8),c(coef(rif.pois)[3],coef(rif.qpois)[3],coef(rif.negbin)[3],coef(rif.zip)[3],coef(rif.hur)[3]),pch=21,bg=c("black","red","blue","orange","brown"),col=c("black","red","blue","orange","brown"))						

# effet SHA500
arrows(x0=c(5,5.2,5.4,5.6,5.8),x1=c(5,5.2,5.4,5.6,5.8),
			 y0=c(conf.pois[4,1],conf.qpois[4,1],conf.negbin[4,1],conf.zip[4,1],conf.hur[4,1])
			 ,y1=c(conf.pois[4,2],conf.qpois[4,2],conf.negbin[4,2],conf.zip[4,2],conf.hur[4,2])
			 ,length=0,code=2,lwd=3,col=c("black","red","blue","orange","brown"))
points(c(5,5.2,5.4,5.6,5.8),c(coef(rif.pois)[4],coef(rif.qpois)[4],coef(rif.negbin)[4],coef(rif.zip)[4],coef(rif.hur)[4]),pch=21,bg=c("black","red","blue","orange","brown"),col=c("black","red","blue","orange","brown"))						


# effet r�gion (contraste vs Banks)
arrows(x0=c(7,8,9),x1=c(7,8,9),
			 y0=conf.zip[6:8,1]
			 ,y1=conf.zip[6:8,2]
			 ,length=0,code=2,lwd=3,col="orange")
points(c(7,8,9),coef(rif.zip)[6:8],col="orange",bg="orange",pch=21)

arrows(x0=c(7.2,8.2,9.2),x1=c(7.2,8.2,9.2),
			 y0=conf.hur[6:8,1]
			 ,y1=conf.hur[6:8,2]
			 ,length=0,code=2,lwd=3,col="brown")
points(c(7.2,8.2,9.2),coef(rif.hur)[6:8],col="brown",bg="brown",pch=21)

legend(x=5,y=4,legend=c("GLM Poisson", "GLM quasi-Poisson","GLM n�g. binom.","ZIP","Hurdle"),lty=rep("solid",5),col=c("black","red","blue","orange","brown"))


################# FIN DU COURS DE MARDI #################################
