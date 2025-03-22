#---------------------------------------------------------------------#
### Formation "certification statistique pour les écologues" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# ce script permet la réplication des exemples du cours du lundi 
#---------------------------------------------------------------------#

# packages nécessaires (certains seront re-chargés par la suite pour la clarté du script mais ce n'est pas nécessaire)
# attention : il y a un conflit entre sjPlot et cowplot pour la fonction plot_grid. Afin d'obtenir le comportement désiré
# de cette fonction dans le cadre de ce script, saisir : cowplot::plot_grid()
library(ggplot2)
library(formattable)
library(cowplot)
library(sjPlot)
library(effects)

# répertoire de travail (à modifier en fonction de votre chemin d'accès)
setwd("C:/Users/jeany/OneDrive/Documents/certificat_2023/niveau2/supports/donnees")

#----------------------------#
### chargement des données ###
#----------------------------#

# extrait de données de Prodon et al, Global Change Biology, 2017
amphi=read.table("amphibiens.txt",header=T,sep="\t")

#----------------------------------------------#
### tendance temporelle sur une seule espèce ###
#----------------------------------------------#

# triton marbré
trimar=subset(amphi,codesp=="TRIMAR")

# exploration données brutes
ggplot(trimar)+
	geom_point(aes(x=an,y=julian))+
	xlab("années")+
	ylab("Date julienne")

# histogramme des dates
hist(trimar$julian,xlab="dates de sortie",ylab="fréquence",main="",cex.lab=1.5,col="gray80")

# tendance
trimar.trend=lm(julian~an,data=trimar)

# résidus
par(mfrow=c(2,2))
plot(trimar.trend)

diff.sp = lm(julian~factor(codesp),data = amphi)
par(mfrow=c(2,2))
plot(diff.sp)
summary(diff.sp)
# résultats
tab_model(trimar.trend,show.se=T,show.stat=T,show.df=T)

amphi$codesp = factor(amphi$codesp,levels=c("BUFBUF","HYLMER","RANTEM","TRIMAR","ALYOBS"))
diff.sp = lm(julian~(codesp),data = amphi)

# sortie graphique
ggplot(trimar)+
	aes(x=an,y=julian)+
		geom_point()+
			geom_smooth(method='lm')+
				xlab("années")+
					ylab("Date julienne")
				
#-------------------------#
### régression multiple ###
#-------------------------#


plot(amphi$an,amphi$alt)

# sur tous les amphibiens de ce jeu de données
ggplot(amphi)+
	geom_point(aes(x=an,y=julian))+
	xlab("années")+
	ylab("Date julienne")

# distribution des dates de sortie
hist(amphi$julian,xlab="dates de sortie",ylab="fréquence",main="",cex.lab=1.5,col="gray80")

# effet espèces
ggplot(amphi)+
	geom_boxplot(aes(x=codesp,y=julian))+
	xlab("espèces")+
	ylab("Date julienne")

# effet altitude
ggplot(amphi)+
	geom_point(aes(x=alt,y=julian))+
	xlab("altitude (m)")+
	ylab("Date julienne")

# régressions linéaires simples
lm.an=lm(julian~an,data=amphi)
tab_model(lm.an,show.se=T,show.stat=T,show.df=T)

lm.alt=lm(julian~alt,data=amphi)
tab_model(lm.alt,show.se=T,show.stat=T,show.df=T)

lm.esp=lm(julian~codesp,data=amphi)
tab_model(lm.esp,show.se=T,show.stat=T,show.df=T)

# sortie graphique
ggplot(amphi)+
	aes(x=an,y=julian)+
		geom_point()+
			geom_smooth(method='lm')+
				xlab("années")+
					ylab("Date julienne")

# régression linéaire multiple
lm.mult=lm(julian~an+alt+codesp,data=amphi)

# résidus
par(mfrow=c(2,2))
plot(lm.mult)

# résumé
summary(lm.mult)

# intervalles de confiance
confint(lm.mult)

# table de résultats
tab_model(lm.mult,show.se=T,show.stat=T,show.df=T)

#----------------------------------------------------------#
### représentations graphiques de la régression multiple ###
#----------------------------------------------------------#

# graphique incorrect pour les années (pour illustration de ce qu'il ne faut pas faire, à ne pas réutiliser)
ggplot(amphi)+
	aes(x=an,y=julian)+
		geom_point()+
			geom_abline(intercept=coef(lm.mult)[1],slope=coef(lm.mult)[2])

# graphique incorrect pour les altitudes (pour illustration de ce qu'il ne faut pas faire, à ne pas réutiliser)
ggplot(amphi)+
	aes(x=alt,y=julian)+
	geom_point()+
	geom_abline(intercept=coef(lm.mult)[1],slope=coef(lm.mult)[3])


# graphique correct avec les résidus partiels (plusieurs solutions parmi de nombreuses disponibles)
library(jtools)
p1=effect_plot(lm.mult, pred = an, interval = TRUE, partial.residuals = TRUE)
p2=effect_plot(lm.mult, pred = alt, interval = TRUE, partial.residuals = TRUE)
p3=effect_plot(lm.mult, pred = codesp, interval = TRUE, partial.residuals = TRUE)
cowplot::plot_grid(p1,p2,p3)

library(ggeffects) # selon la fonction utilisée, le comportement de ggeffects varie
p1=ggpredict(lm.mult,terms="an") # ggpredict : effet de la variable conditionné aux valeurs de référence (codesp=ALYOBS et alt=0)
p1.alt=ggemmeans(lm.mult,terms="an") # ggemmeans : effet marginal à altitude moyenne et pour une espèce moyenne
p1b=plot(p1,residuals=T)
p1.altb=plot(p1.alt,residuals=T)
cowplot::plot_grid(p1b,p1.altb)

p2=ggpredict(lm.mult,terms="alt")
p3=ggpredict(lm.mult,terms="codesp")
p2b=plot(p2,residuals=T)
p3b=plot(p3,residuals=T)
cowplot::plot_grid(p1b,p2b,p3b)

library(sjPlot)
plot_model(lm.mult,terms="an",type="pred")
plot_model(lm.mult,terms="an",type="pred",show.data=T) # attention, représente les données brutes et non les résidus partiels (donc incorrect dans la majorité des cas)!
plot_model(lm.mult,terms="alt",type="pred",show.data=T) # attention, représente les données brutes et non les résidus partiels (donc incorrect dans la majorité des cas)! 

library(visreg) # solution en graphiques R de base, pratique car très flexible (de nombreux types de modèles implémentés dont modèles mixtes, régressions quantiles, etc)
par(mfrow=c(2,2))
visreg(lm.mult,xvar="an")
visreg(lm.mult,xvar="alt")
visreg(lm.mult,xvar="codesp")

library(effects) # autre solution en graphiques R de base flexible
eff.mod=effect(term="an",mod=lm.mult)
plot(eff.mod) # sans les résidus partiels
eff.mod2=predictorEffects(lm.mult, residuals=TRUE)
plot(eff.mod2,partial.residuals=list(smooth=F)) # avec les résidus partiels

#----------------------------------------------------#
### le même modèle sur variables centrées-réduites ###
#----------------------------------------------------#

# la fonction scale() permet de centrer réduire (ou de juste centrer, ou juste réduire, selon les arguments choisis)
# à noter qu'on ne centre-réduit pas la variable codesp qui est un facteur. C'est néanmoins possible mais plus complexe
# et pas forcément intéressant pour notre objectif ici. 
lm.mult.sc=lm(julian~scale(an)+scale(alt)+codesp,data=amphi)
summary(lm.mult.sc)
tab_model(lm.mult.sc,show.se=T,show.stat=T,show.df=T)

# représentation utile pour comparer des effets de variables centrées réduites, avec sjPlot
plot_model(lm.mult.sc,terms=c("scale(an)","scale(alt)"),type="est")

###### FIN DU SCRIPT DU COURS DE LUNDI ########