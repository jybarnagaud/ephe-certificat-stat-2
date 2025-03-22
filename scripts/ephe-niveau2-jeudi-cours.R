#---------------------------------------------------------------------#
### Formation "certification statistique pour les écologues" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# ce script permet la réplication des exemples du cours du mercredi 
#---------------------------------------------------------------------#
library(ggplot2)
library(ggeffects)
library(visreg)
library(interactions)
library(cowplot)
library(investr)

setwd("D:/certificat_2023/niveau2/supports/donnees")

#------------------------------------------------------------------#
### Théorie des îles ou théorie de la niche? Comparer des modèles###
#------------------------------------------------------------------#

# test des ratios de vraisemblance (LRT)

# charger les données
lomolino=read.table("Lomolino.txt",header=T,sep="\t")
lomolino$log.SURF=log(lomolino$SURF) # on log transforme la surface, voir cours sur le GLM pour justification

# 2 hypothèses emboîtées
m.spatial=glm(NSPEC~log.SURF+LATI+DSOUR,family="poisson",data=lomolino)
m.surface=glm(NSPEC~log.SURF,family="poisson",data=lomolino)

# on contrôle les résidus du premier modèle
par(mfrow=c(2,2))
plot(m.spatial)

# et du deuxième
par(mfrow=c(2,2))
plot(m.surface)

# vraisemblances (log transformées)
logLik(m.spatial)
logLik(m.surface)

# ratio de vraisemblances
library(lmtest)
lrtest(m.spatial,m.surface)

## Comparer des modèles non emboîtés

# 2 hypothèses : une sur des variables spatiales, une sur des variables écologiques
m.spatial=glm(NSPEC~log.SURF+LATI+DSOUR,family="poisson",data=lomolino)
m.ecol=glm(NSPEC~DFOR+DGRASS+DPLEIS+LATI,family="poisson",data=lomolino)

# on contrôle les résidus
par(mfrow=c(2,2))
plot(m.spatial)

# l'autre modèle
par(mfrow=c(2,2))
plot(m.ecol)

# vraisemblances
logLik(m.spatial)
logLik(m.ecol)

# AIC
AIC(m.spatial)
AIC(m.ecol)

# exemple fictif sur l'AIC : explication du critère de 2 unités d'AIC
L1=10 # vraisemblance du modèle 1
L2=10 # vraisemblance du modèle 2

AIC1=2*3-2*log(10)
AIC2=2*4-2*log(10)

# on compare les deux hypothèses (données de Lomolino) avec l'AIC
AIC(m.spatial,m.ecol)

summary(m.spatial) # on a déjà étudié ce modèle dans le cours sur le GLM

#------------------------------------------------#
### Phénologie comparée des chênes et des pins ###
# une interaction entre une variable catégorique
# et une variable continue
#------------------------------------------------#

# les données : bourgeonnement pour 2 espèces : pin sylvestre et chêne sessile (scores de 1 = non débourré à 9 = feuille complète)
brg=read.table("phenologie_arbres.txt",header=T,sep="\t")

# exploration des données
ggplot(brg)+
	aes(x=date_jul,y=bmean,color=essence)+
		geom_point(size=2)+
			labs(title="",x="date julienne",y="degré de débourrement")+
				scale_color_manual(values = c("#E69F00", "#0072B2"),labels=c("Chêne","Pin"))
						

# deux modèles séparés
bg.ch=lm(bmean~date_jul,data=brg,subset=essence=="CHS")
bg.ps=lm(bmean~date_jul,data=brg,subset=essence=="PS")

ggplot(brg)+
	aes(x=date_jul,y=bmean,color=essence)+
	geom_point(size=2)+
	labs(title="",x="date julienne",y="degré de débourrement")+
		scale_color_manual(values = c("#E69F00", "#0072B2"),labels=c("Chêne","Pin"))+
			geom_abline(intercept=coef(bg.ch)[1],slope=coef(bg.ch)[2],color="#E69F00")+
			geom_abline(intercept=coef(bg.ps)[1],slope=coef(bg.ps)[2],color="#0072B2")
	
				

# on a l'impression que les pentes différent mais on n'en est pas sûr --> on teste l'interaction plutôt que séparer 2 modèles
# modèles (on montre les deux structures)
bg.mod=lm(bmean~date_jul+essence,data=brg)
bg.mod.int=lm(bmean~date_jul+essence+(date_jul:essence),data=brg)
bg.mod.int=lm(bmean~date_jul*essence,data=brg) 

# résidus
par(mfrow=c(2,2))
plot(bg.mod.int) # on note que c'est pas très propre mais il n'y a pas d'hétéroscédasticité

# on choisit la structure
AIC(bg.mod,bg.mod.int)

# résumé du modèle
summary(bg.mod.int)

# effets marginaux
p1=ggpredict(bg.mod.int,terms=c("date_jul","essence"))
plot(p1,residuals=F,add.data=T)+
	labs(title="",x="date julienne",y="degré de débourrement")+
		scale_color_manual(values = c("#E69F00", "#0072B2"),labels=c("Chêne","Pin"))

# version visreg
visreg(bg.mod.int,xvar="date_jul",by="essence")

# la construction de la représentation du modèle à la main (en utilisant les graphiques R de base)
date.pred=seq(from=min(brg$date_jul),to=max(brg$date_jul),length.out=100)
xpred.ch=list(date_jul=date.pred,essence=rep("CHS",100))
ypred.ch=predict(bg.mod.int,newdata=xpred.ch,interval="confidence")
dev.off()
plot(brg$date_jul,brg$bmean,col=as.numeric(brg$essence),xlab="Date Julienne",ylab="Bourgeonnement",bg=as.numeric(brg$essence),pch=21)
lines(date.pred,ypred.ch[,1])
lines(date.pred,ypred.ch[,2],lty="dashed")
lines(date.pred,ypred.ch[,3],lty="dashed")

xpred.ps=list(date_jul=date.pred,essence=rep("PS",100))
ypred.ps=predict(bg.mod.int,newdata=xpred.ps,interval="confidence")
lines(date.pred,ypred.ps[,1],col="red")
lines(date.pred,ypred.ps[,2],lty="dashed",col="red")
lines(date.pred,ypred.ps[,3],lty="dashed",col="red")
legend("bottomright",lty=c("solid","solid"),col=c("black","red"),legend=c("Chêne","Pin"))

#--------------------------------------------------------#
### Effet du paysage sur la chenille processionnaire : ###
# une interaction entre deux variables continues
#--------------------------------------------------------#

# données de chenilles déjà exploitées pendant le cours sur le GLM
chenilles=read.table("chenilles.txt",header=T,sep="\t")

# proportion d'attaques
chenilles$prop_attaq=chenilles$nbattaq/chenilles$nbpins

# effets joints de prop_pin (proportion de pins dans la forêt) et de prop_foret (proportion de forêt dans le paysage)
col.prop=brewer.pal(10,"Blues")

ggplot(chenilles)+
	aes(x=prop_pin,y=prop_foret,fill=prop_attaq)+
		geom_point(size=4,pch=21,alpha=0.8)+
				scale_fill_gradient(low="#F7FBFF",high="#08306B")+
					labs(title="",x="Proportion de pins par parcelle, %",y="Proportion de forêt dans le paysage, %",fill="Taux d'infestation, %")
			
ggplot(chenilles)+
	aes(x=prop_pin,y=prop_attaq,fill=prop_foret)+
	geom_point(size=4,pch=21,alpha=0.8)+
	scale_fill_gradient(low="seashell",high="green4")+
	labs(title="",x="Proportion de pins par parcelle, %",y="Taux d'infestation, %",fill="Proportion de forêt \n dans le paysage, %")

# modèles : avec et sans interactions
mod.chen=glm(cbind(nbattaq,nbpins-nbattaq)~prop_pin+prop_foret,family="binomial",data=chenilles)
mod.chen.int=glm(cbind(nbattaq,nbpins-nbattaq)~prop_pin*prop_foret,family=binomial,data=chenilles)

# on vérifie les résidus systématiquement
par(mfrow=c(2,2))
plot(mod.chen.int)

# on compare
AIC(mod.chen,mod.chen.int)

# résumé du modèle
summary(mod.chen.int)

# pour la représentation visuelle il va falloir scinder la variable "prop foret" en plusieurs, sjPlot est pratique pour ça
# on peut choisir quelles valeurs on utilise pour couper la variable en changeant mdrt.values, !! ça change l'impression de résultat 

library(sjPlot)
plot_model(mod.chen.int, type = "int",mdrt.values="meansd") # moyenne et sd de la variable
plot_model(mod.chen.int, type = "int",mdrt.values="quart") # quartiles
plot_model(mod.chen.int, type = "int",mdrt.values="minmax") # min et max (à éviter, ne donnent que les extrêmes, il y a des méthodes stats plus appropriées

## Dynamiques prairiales entre sites : représenter une interactions entre deux facteurs
#--------------------------------------------------------------------------------------

calais=read.table("vegetation_pas_de_calais.txt",header=T,sep="\t")
calais$time.period=relevel(calais$time.period,ref="fin_90")

ggplot(calais)+
		aes(x=time.period,y=richesse)+
				geom_boxplot()+
					facet_wrap(~Secteur,ncol=2)

flore.calais0=lm(hveg~time.period+Secteur,data=calais)
par(mfrow=c(2,2))
plot(flore.calais0)
summary(flore.calais0)

flore.calais=lm(hveg~time.period*Secteur,data=calais)

par(mfrow=c(2,2))
plot(flore.calais)
summary(flore.calais) 
# pas de variation différenciée par secteur malgré des différences entre secteurs = ces différences sont restées les mm entre périodes

# représentation graphique (une parmi d'autres)
interactions::cat_plot(flore.calais,pred=Secteur,modx=time.period,interval=T,colors=c("blue4","tan3"),interval.geom="linerange")+
	labs(title="",x="Secteur",y="Hauteur végétation estimée, cm")

#--------------------------------------------------------------------------#
### 9.3 Modéliser un optimum de réponse le long d’un gradient écologique ###
#--------------------------------------------------------------------------#

## L’effet du climat sur la reproduction des tétras : un modèle à effets quadratiques

# données tétras préalpes
tetras=read.table("tetras_lyre.txt",header=T,sep="\t")
tetras$prop=tetras$Nichees/tetras$Poules
tetras$UN=factor(tetras$UN)

summary(tetras)

# exploration graphique

# NAO
synth.nao=cbind(unique(tetras$Annee),
tapply(tetras$NAOdjfm,INDEX=tetras$Annee,FUN="mean"))
colnames(synth.nao)=c("annee","moy")
synth.nao=as.data.frame(synth.nao)
synth.nao=synth.nao[order(synth.nao$annee),]
plot(synth.nao$annee,synth.nao$moy,type="b",cex=1.5,pch=21,bg="black",xlab="Années",ylab="NA0 (décembre - mars")
abline(h=0,lty="dashed")
		
# hypothèse
library(boot)
nao.sim=rnorm(1000,mean(synth.nao[,2]),sd(synth.nao[,2]))
y=0.5+0.1*nao.sim-0.3*(nao.sim^2)
psi=inv.logit(y)
z=cbind(nao.sim,psi)
z=z[order(z[,1]),]
plot(z[,1],z[,2],type="l",lwd=2,xlab="NAO décembre - mars",ylab="Taux de reproduction des tétras",cex.lab=1.5)
abline(v=0,lty="dashed")

# données
ggplot(tetras)+
	aes(x=NAOdjfm,y=prop)+
		geom_point(size=3)+
			labs(x="NAO décembre - mars",y="Taux de reproduction des tétras")

# modèles
tetras.mod.lin=glm(cbind(Nichees,Poules-Nichees)~NAOdjfm+RN,family="binomial",data=tetras)
tetras.mod.quad=glm(cbind(Nichees,Poules-Nichees)~poly(NAOdjfm,2,raw=T)+RN,family="binomial",data=tetras)
tetras.mod.quad=glm(cbind(Nichees,Poules-Nichees)~NAOdjfm+I(NAOdjfm^2)+RN,family="binomial",data=tetras)
tetras.mod.quad=glm(cbind(Nichees,Poules-Nichees)~poly(NAOdjfm,2)+RN,family="binomial",data=tetras)

AIC(tetras.mod.lin,tetras.mod.quad)

# résidus
par(mfrow=c(2,2))
plot(tetras.mod.lin)

par(mfrow=c(2,2))
plot(tetras.mod.quad)

# sortie du modèle
summary(tetras.mod.quad)

# représentation des effets
p1=plot(ggpredict(tetras.mod.quad,terms=c("NAOdjfm")),residuals=T)
p2=plot(ggpredict(tetras.mod.quad,terms=c("RN")),residuals=T)
plot_grid(p1,p2)


# trouver l'optimum de réponse 

y=function(x){coef(tetras.mod.quad)[1]+coef(tetras.mod.quad)[2]*x+coef(tetras.mod.quad)[3]*x*x}
xmax=optimize(y, interval=range(tetras$NAOdjfm), maximum=TRUE)
xmax

# autre manière de faire
library(polynom)
p0 <- polynomial(c(coef(tetras.mod.quad)[1],coef(tetras.mod.quad)[2],coef(tetras.mod.quad)[3]))
p1 <- deriv(p0); p2 <- deriv(p1) # dérivées 1 et 2   
xm <- solve(p1) # max et min de p0
xmax = xm[predict(p2, xm) < 0] # select the maxima 
xmax 

# plot
p1=ggpredict(tetras.mod.quad,terms=c("NAOdjfm"))
plot(p1,residuals=T)+
	geom_vline(xintercept=xmax,linetype="dotted",colour="steelblue",size=2)+
		labs(title="",y="% nichées",x="NAO décembre - mars")


## terme d'interaction
ggplot(tetras)+
	aes(x=NAOdjfm,y=prop)+
		geom_point()+
			facet_wrap(~RN,nrow=2)

# faire avec un terme d'interaction - UN
tetras.mod.quad.int=glm(cbind(Nichees,Poules-Nichees)~poly(NAOdjfm,2,raw=T)*RN,family="binomial",data=tetras)
par(mfrow=c(2,2))
plot(tetras.mod.quad.int)

# représentation graphique
tetras.mod.quad.int=glm(cbind(Nichees,Poules-Nichees)~(NAOdjfm+I(NAOdjfm^2))*RN,family="binomial",data=tetras)
p1=ggpredict(tetras.mod.quad.int,terms=c("NAOdjfm","RN"))
plot(p1,residuals=T,facets=T,colors="bw")+
	labs(title="",x="NAO décembre - mars",y="%Nichées")

# attention néanmoins à ne pas abuser des interactions : 
AIC(tetras.mod.quad.int,tetras.mod.quad)
# le modèle est manifestement surfitté
summary(tetras.mod.quad.int)


#-----------------------------------------------------#
### Représenter une variation non linéaire complexe ###
#-----------------------------------------------------#

## Une variation temporelle complexe

amphi=read.table("amphibiens.txt",header=T,sep="\t")

# d'abord un modèle quadratique
rantem=subset(amphi,codesp=="RANTEM")
ggplot(rantem)+
	aes(x=an,y=julian)+
		geom_point()+
			labs(title="",x="Années",y="Date de sortie d'hibernation")

# modèle linéaire
rantem.mod=lm(julian~an,data=rantem)
par(mfrow=c(2,2))
plot(rantem.mod)

# avec effet quadratique
rantem.mod.quad=lm(julian~poly(an,2,raw=T),data=rantem)
par(mfrow=c(2,2))
plot(rantem.mod.quad)

AIC(rantem.mod,rantem.mod.quad)

# représentation
rantem.mod.quad=lm(julian~an+I(an^2),data=rantem)
p1=ggemmeans(rantem.mod.quad,terms="an")
plot(p1,add.data=T)+
	labs(title="",x="Années",y="Date de sortie d'hibernation")

# avec un GAM
library(mgcv)
rantem.mod.gam=gam(julian~s(an),data=rantem)
gam.check(rantem.mod.gam)
abline(0,1,col="red")

# vérification des résidus et de l'ajustement
gam.check(rantem.mod.gam,old.style=T)
plot(predict(rantem.mod.gam),residuals(rantem.mod.gam),xlab="valeurs prédites",ylab="résidus",main="Residuals vs fitted",font.main=2)
abline(h=0,lty="dashed")

summary(rantem.mod.gam)

plot(rantem.mod.gam)
AIC(rantem.mod.gam,rantem.mod.quad,rantem.mod)

# comparaison graphique
dev.off()
plot(rantem.mod.gam,shift=coef(rantem.mod.gam)[1],se=F,ylim=range(rantem$julian),lwd=2,rug=F,xlab="années",ylab="date julienne",main="Grenouille rousse")
points(rantem$an,rantem$julian,pch=21,bg="gray70",col="gray70",cex=1.5)
p.quad=predict(rantem.mod.quad,newdata=list(an=seq(min(rantem$an),max(rantem$an),by=1)))
lines(seq(min(rantem$an),max(rantem$an),by=1),p.quad,col="red",lwd=2)
legend("topright",lwd=2,col=c("black","red"),legend=c("GAM","modèle quadratique"))


# on contraint les k pour montrer l'effet du nombre de noeuds
rantem.mod.gam.k1=gam(julian~s(an,k=1),data=rantem)
rantem.mod.gam.k2=gam(julian~s(an,k=2),data=rantem)
rantem.mod.gam.k3=gam(julian~s(an,k=3),data=rantem)
rantem.mod.gam.k4=gam(julian~s(an,k=4),data=rantem)
rantem.mod.gam.k5=gam(julian~s(an,k=5),data=rantem)
rantem.mod.gam.k6=gam(julian~s(an,k=6),data=rantem)

par(mfrow=c(2,2))
#plot(rantem.mod.gam.k1)
#plot(rantem.mod.gam.k2)
plot(rantem.mod.gam.k3,residuals=T,cex=5,rug=F,shade=T,ylab="Tendance estimée")
legend("topright",legend=paste("edf =",round(summary(rantem.mod.gam.k3)$edf,2),sep=" "),bty="n")
plot(rantem.mod.gam.k4,residuals=T,cex=5,rug=F,shade=T,ylab="Tendance estimée")
legend("topright",legend=paste("edf =",round(summary(rantem.mod.gam.k4)$edf,2),sep=" "),bty="n")
plot(rantem.mod.gam.k5,residuals=T,cex=5,rug=F,shade=T,ylab="Tendance estimée")
legend("topright",legend=paste("edf =",round(summary(rantem.mod.gam.k5)$edf,2),sep=" "),bty="n")
plot(rantem.mod.gam.k6,residuals=T,cex=5,rug=F,shade=T,ylab="Tendance estimée")
legend("topright",legend=paste("edf =",round(summary(rantem.mod.gam.k6)$edf,2),sep=" "),bty="n")

# plusieurs splines dans un modèle
multi.gam=gam(julian~s(an)+s(alt),data=rantem)
summary(multi.gam)
par(mfrow=c(1,2))
plot(multi.gam)

# une spline et un effet linéaire
multi.gam.lin=gam(julian~s(an)+alt,data=rantem)
AIC(multi.gam,multi.gam.lin)

summary(multi.gam)


# et avec des effets linéaires en plus
multi.gam2=gam(julian~s(an)+s(alt)+codesp,data=amphi)
summary(multi.gam2)
par(mfrow=c(2,2))
visreg(multi.gam2)

# il vaut peut être mieux contraindre un peu s(an) qui semble surlisser : est-ce qu'on perd vraiment en imposant k=3?
multi.gam3=gam(julian~s(an,k=3)+s(alt,k=3)+codesp,data=amphi)
summary(multi.gam2)
par(mfrow=c(2,2))
visreg(multi.gam3)

AIC(multi.gam2,multi.gam3) # numériquement oui, mais en termes d'interprétation ça ne change vraiment rien

# un gam toutes espèces
amphi$codesp = factor(amphi$codesp)

all.sp.gam=gam(julian~codesp+s(an)+s(alt)+s(an,by=codesp)+s(alt,by=codesp),data=amphi)
AIC(all.sp.gam,multi.gam2)

gam.check(all.sp.gam)
summary(all.sp.gam)

par(mfrow=c(2,5))
plot(all.sp.gam,residuals=T)

#------------------------------------------------------------------------------------#
### Représenter la phénologie du bourgeonnement plus finement : le modèle sigmoïde ###
#------------------------------------------------------------------------------------#

deb.brg0=read.table("phenologie_arbres.txt",header=T,sep="\t")
deb.brg=subset(deb.brg0,!is.na(bmean))

deb.brg.pin=subset(deb.brg,essence=="PS" )
deb.brg.ch=subset(deb.brg,essence=="CHS" )

# données brutes
ggplot(deb.brg)+
	aes(x=factor(date_jul),y=bmean)+
		geom_boxplot(aes(fill=essence))+
			scale_fill_manual(values=c("darkgreen","darkred"),labels=c("Chêne","Pin"))+
				labs(x="Date julienne",y="Débourrement moyen par arbres",fill="Essence")

				
# simul nls
xx=rnorm(1000,2,1)
y=function(x,a,b,c,d){
	a+((b-a)/(1+exp((c-x)/d)))
	}

z=cbind(xx,y(xx,1,5,2,0.5))
z=z[order(z[,1]),]
plot(z,type="l",xlab="X",ylab="sigmoïde(x)",lwd=2)

# d'abord les pins
deb.pin.mod=nls(bmean~SSfpl(date_jul,a,b,c,d),data=deb.brg.pin)
summary(deb.pin.mod)

# puis les chênes
deb.ch.mod=nls(bmean~SSfpl(date_jul,a,b,c,d),data=deb.brg.ch)
summary(deb.ch.mod)

# Représentation graphique avec le package investr (plus compliqué sous ggplot)
plotFit(deb.pin.mod, interval = "prediction", pch = 19
				, col.pred = adjustcolor("darkorange4", 0.5),xlim=c(105,155), ylim=c(0,10),shade = TRUE,xlab="Date julienne",ylab="Degré de bourgeonnement")

par(new = TRUE)
plotFit(deb.ch.mod, interval = "prediction", pch = 19, 
				col.pred =  adjustcolor("chartreuse4", 0.5),xlim=c(105,155), ylim=c(0,10), shade = TRUE,xlab="Date julienne",ylab="Degré de bourgeonnement")

legend("bottomright",legend=c("Pin sylvestre","Chêne sessile"),fill=c("darkorange4","chartreuse4"),bty="n")

# le même graphique à la main pour illustrer la construction des intervalles de confiance (fonction predict)
prd.pin=predict(deb.pin.mod)
z=cbind(prd.pin,deb.brg.pin$date_jul,as.numeric(factor(deb.brg.pin$date_jul)))
z=z[order(z[,3]),]
boxplot(deb.brg.pin$bmean~deb.brg.pin$date_jul,col="brown",ylim=range(deb.brg$bmean),xlab="Date Julienne",ylab="Degré de débourrement")
boxplot(deb.brg.ch$bmean~deb.brg.ch$date_jul,col="darkgreen",add=T)
legend("bottomright",legend=c("Pin sylvestre","Chêne sessile"),fill=c("brown","darkgreen"))
lines(z[,3],z[,1],lwd=2,col="red")

prd.ch=predict(deb.ch.mod)
z=cbind(prd.ch,deb.brg.ch$date_jul,as.numeric(factor(deb.brg.ch$date_jul)))
z=z[order(z[,3]),]
lines(z[,3],z[,1],lwd=2,col="darkblue")

# prédire les valeurs manquantes dans deb.brg0
brg.na.ch=subset(deb.brg0,is.na(bmean) & essence=="CHS")
pred.ch=predFit(deb.ch.mod,newdata=data.frame(date_jul=brg.na.ch$date_jul),interval="prediction")
# noter le choix de l'intervalle de confiance "prediction" (plutôt que "confidence" qui sert à estimer l'IC de la courbe)
# on utilise la fonction predFit de investr plutôt que predict.nls, qui ne calcule pas les intervalles de confiance (les arguments 
# correspondants sont ignorés). La fonction predFit fait une approximation linéaire (= delta method) pour calculer les IC

brg.na.pin=subset(deb.brg0,is.na(bmean) & essence=="PS")
pred.pin=predFit(deb.pin.mod,newdata=data.frame(date_jul=brg.na.pin$date_jul),interval="prediction")

# ajout des valeurs prédites à la représentation graphique
plotFit(deb.pin.mod, interval = "prediction", pch = 19
				, col.pred = adjustcolor("darkorange4", 0.5),xlim=c(105,155), ylim=c(0,10),shade = TRUE,xlab="Date julienne",ylab="Degré de bourgeonnement")
par(new = TRUE)
plotFit(deb.ch.mod, interval = "prediction", pch = 19, 
				col.pred =  adjustcolor("chartreuse4", 0.5),xlim=c(105,155), ylim=c(0,10), shade = TRUE,xlab="Date julienne",ylab="Degré de bourgeonnement")
legend("bottomright",legend=c("Pin sylvestre","Chêne sessile"),fill=c("darkorange4","chartreuse4"),bty="n")

points(brg.na.ch$date_jul,pred.ch[,"fit"],pch=3,col="darkgreen",cex=1.5)
segments(x0=brg.na.ch$date_jul,x1=brg.na.ch$date_jul,y0=pred.ch[,"lwr"],y1=pred.ch[,"upr"],col="darkgreen")

points(brg.na.pin$date_jul,pred.pin[,"fit"],pch=3,col="red4",cex=1.5)
segments(x0=brg.na.pin$date_jul,x1=brg.na.pin$date_jul,y0=pred.pin[,"lwr"],y1=pred.pin[,"upr"],col="red4")
