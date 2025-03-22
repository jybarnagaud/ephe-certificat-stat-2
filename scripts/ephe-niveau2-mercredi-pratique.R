#---------------------------------------------------------------------------#
### Formation "Analyse de données pour les écologues renforcement" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# TD du mercredi : modèle linéaire généralisé
#---------------------------------------------------------------------------#

library(ggplot2)
library(arm)
library(questionr)
library(ggeffects)
library(MASS)
library(corrplot)

rm(list=ls())
setwd("D:/certificat_2023/niveau2/supports/donnees")

#----------------------------------------------------------------#
### Cas pratique 1 : Ecologie comportementale sur les baleines ###
#----------------------------------------------------------------#

# charger les données
whale=read.table("baleines.txt",header=T,sep="\t")

# on s'intéresse au temps passé sur le stimulus
hist(whale$time_z12)

# on ne peut pas garder cette variable telle quelle : en première intention, on pourrait vouloir
# se contenter d'une transformation en log ou en racine carrée, mais si on réfléchit un peu à ce qu'elle
# veut dire, il y a clairement 2 process séparés : d'une part, est-ce que le groupe vient vers le stimulus
# ou est-ce qu'il se contente de transiter, d'autre part combien de temps il passe près du stimulus lorsqu'il
# est venu. 

# on représente le fait que le groupe se dirige ou non vers le stimulus par une variable binaire
whale$is.near=NA
whale[which(whale$time_z12>0),"is.near"]=1 # le groupe vient vers le stimulus s'il passe au moins >0 min près du dispositif
whale[which(whale$time_z12==0),"is.near"]=0 # le groupe ne vient pas vers le stimulus s'il passe 0 min près du dispositif

# si je viens, combien de temps je reste : cette question n'a de sens que si le temps passé est >0
whale2=subset(whale,time_z12>0)
hist(whale2$time_z12)

#--------------------------------------------------------------------------------------#
### analyse de la variable "is.near" : les baleines viennent ou non vers le stimulus ###
#--------------------------------------------------------------------------------------#


## Hypothèse 1 : effet attractif du Krill et du DMS, mais pas des deux témoins
## Hypothèse 2 : effet attractif plus fort dans les zones d'alimentation que de reproduction
## Hypothèse 3 : effet d'entraînement : les gros groupes répondent plus


# nombre de fois qu'un groupe s'intéresse / ne s'intéresse pas au dispositif pour chaque stimulus
is.near.stim=as.data.frame(table(whale[,c("is.near","chemical")]))
ggplot(data=is.near.stim, aes(x=chemical, y=Freq, fill=is.near)) + geom_bar(stat="identity",position="dodge")

# nombre de fois qu'un groupe s'intéresse / ne s'intéresse pas au dispositif dans chaque zone
is.near.pl=as.data.frame(table(whale[,c("is.near","place")]))
ggplot(data=is.near.pl, aes(x=place, y=Freq, fill=is.near)) + geom_bar(stat="identity",position="dodge")

# effet taille de groupe
ggplot(whale)+aes(x=factor(is.near),y=nb_indiv)+geom_boxplot() # il semble y avoir un effet d'entraînement, les individus seuls ne s'intéressent pas

# avant de modéliser, on va changer l'ordre des facteurs de chemical afin que le contrôle soit le niveau de référence
whale$chemical=factor(whale$chemical,levels=c("CTL","CLAY","DMS","KRILL"))

# et on met la seule zone de reproduction comme niveau de référence
whale$place=factor(whale$place,levels=c("Madagascar","Iceland","Antarctica"))

# GLM binomial : modélise la probabilité de s'intéresser au dispositif en fonction du stimulus, de l'endroit et du nombre d'individus
glm.baleines=glm(is.near~chemical+place+nb_indiv,family="binomial",data=whale)

# notez que comme la variable est binaire, il n'est pas utile de spécifier le nombre d'essais

# vérification des résidus
par(mfrow=c(2,2))
plot(glm.baleines) # ces résidus sont certes structurés, mais c'est attendu : on voit persister la structure binomiale dans les résidus

# pour vous convaincre que cette structure résiduelle ne pose pas de problème, voici une simulation qui se place dans les conditions idéales
# du GLM : des résidus normaux centrés sur 0, une variable de réponse tirée dans une binomiale. Augmentez la variance de eps pour faire baisser
# le pourcentage de déviance expliquée et observez les résidus. 
a=1
b=0.05
eps=rnorm(100,0,1)
x=rnorm(100,3,2)
ip=a+b*x+eps
p=invlogit(ip)
y=rbinom(100,1,p)
glm.test=glm(y~x,family=binomial)
par(mfrow=c(2,2))
plot(glm.test)

# pourcentage de déviance expliquée
100*(155-122)/155

# inférence sur le modèle
summary(glm.baleines)

# calcul des odds ratios
odds.ratio(glm.baleines)

# on voit que le seul traitement qui a un petit effet est le krill (probabilité de s'approcher multipliée par 3 en présence de krill par 
# rapport au témoin, mais il y a une forte incertitude qui, liée à un faible échantillon, va nous conduire à ne pas conclure sur cet effet.

# par ailleurs, on voit que les baleines s'approchent moins du dispositif dans les deux zones d'alimentation que dans la zone de reproduction
# prise comme référence, ce qui va contre notre prédiction

# enfin, il y a un effet d'entraînement lié au groupe (probabilité de s'approcher multipliée par 3 à chaque individu en plus dans le groupe)

# représentation des effets
p1=plot(ggemmeans(glm.baleines,terms="chemical"),residuals=T)+labs(title="",x="Stimulus",y="Probabilité d'approche marginale")
p2=plot(ggemmeans(glm.baleines,terms="place"),residuals=T)+labs(title="",x="Site",y="Probabilité d'approche marginale")
p3=plot(ggemmeans(glm.baleines,terms="nb_indiv"),residuals=T)+labs(title="",x="nombre d'individus",y="Probabilité d'approche marginale")
cowplot::plot_grid(p1,p2,p3)

#---------------------------------------------------------------------------------------#
### analyse du temps passé près du stimulus sachant que les baleines s'en sont approchées ###
#---------------------------------------------------------------------------------------#

# ensuite, on analyse le temps passé, sachant qu'un groupe est venu
# il s'agit d'une variable quantitative - la distribution est asymétrique mais 
# il y a beaucoup de variabilité sur tout le gradient de temps, donc il n'est peut-être
# pas nécessaire de transformer les données
lm.temps=lm(time_z12~chemical+place+nb_indiv,data=whale2)

# résidus : très léger patron d'hétéroscédasticité, mais le nuage de points ne semble pas si structuré
par(mfrow=c(2,2))
plot(lm.temps)

# inférence
summary(lm.temps)

p1=plot(ggemmeans(lm.temps,terms="chemical"),residuals=T)+labs(title="",x="Stimulus",y="Temps passé marginal")
p2=plot(ggemmeans(lm.temps,terms="place"),residuals=T)+labs(title="",x="Site",y="Temps passé marginal")
p3=plot(ggemmeans(lm.temps,terms="nb_indiv"),residuals=T)+labs(title="",x="Nombre d'individus",y="Temps passé marginal")
cowplot::plot_grid(p1,p2,p3)

# on constate un assez fort effet du krill, qui augmente de 7 minutes le temps passé près du stimulus par rapport au contrôle
# en revanche, il n'y a pas d'effet du nombre d'individus ni du site.

#----------------#
### Conclusion ###
#----------------#

# le fait de s'intéresser ou non au dispositif semble ne dépendre que du groupe (effet d'entrainement, curiosité?)
# en revanche, il y a une réponse assez nette au stimulus une fois près du dispositif
# Est-ce qu'il faut pour autant en conclure qu'il y a une réponse nette à un stimulus olfactif évoquant l'alimentation?
# Les résultats sont compatibles avec cette hypothèse, mais ils sont un peu faible pour affirmer une causalité
# Il faudrait donc répliquer l'expérience afin d'étudier la reproductibilité de cet effet du krill, et coupler
# ces résultats avec une analyse des structures anatomiques susceptibles de servir de récepteurs olfactifs

# Dans l'article original, le problème a été résolu un peu différemment, car les auteurs ne souhaitaient pas analyser séparément
# la probabilité d'entrer dans la zone, et le temps passé près du dispositif. Ils ont donc utilisé une régression tobit, qui est
# un mélange entre une distribution binomiale et une distribution normale, et dont le but est de traiter les excès de 0 d'une variable
# quantitative continue - c'est ni plus ni moins ce que nous avons fait ici, mais nous avons décomposé les étapes.

#-----------------------------------------------------------------------#
### Cas pratique 2 : déterminants de l'abondance des mouettes rieuses ###
#-----------------------------------------------------------------------#

mouette=read.table("mouette_rieuse_hyeres.txt",header=T,sep="\t")

# Exploration graphique des variables
hist(mouette$comptage) # on voit d'emblée qu'il va y avoir une forte surdispersion dans nos comptages

# corrélation des variables explicatives
xx=cor(mouette[,c("sel","temperature","o2","niveau")])
corrplot(xx) # l'oxygène dissous est assez fortement corrélé à la température, elle-même négativement corrélée à la salinité

# on supprime la température de la liste des variables explicatives 
xx=cor(mouette[,c("sel","o2","niveau")])
corrplot(xx) # cette fois les corrélations restent modérées, on devrait avoir suffisamment de variation pour estimer séparément les effets

# variations temporelles des variables explicatives
p1=ggplot(mouette)+aes(x=factor(date),y=sel)+geom_boxplot()
p2=ggplot(mouette)+aes(x=factor(date),y=o2)+geom_boxplot()
p3=ggplot(mouette)+aes(x=factor(date),y=niveau)+geom_boxplot()
plot_grid(p1,p2,p3) # on voit bien les variations saisonnières des variables explicatives - elles ne nous intéressent pas spécialement ici, mais c'est mieux d'en avoir conscience

# GLM de Poisson
mod.mouet=glm(comptage~sel+o2+niveau,data=mouette,family=poisson)

# contrôle des résidus
par(mfrow=c(2,2))
plot(mod.mouet) # on voit bien la contrainte imposée par les 0, qui implique une forte déviation à la normalité - elle n'est pas gênante en soi, mais
# elle est symptômatique d'un problème dans le modèle. Si on va plus loin, on s'aperçoit que tous les résidus qui dévient de la normalité 
# correspondent à des comptages non nuls : cela suggère que le modèle modélise mal une partie de ce qui nous intéresse 
col=mouette$comptage
col[which(col>0)]=1
col=factor(col)
plot(qqnorm(residuals(mod.mouet)),col=col)

# surdispersion
res.mouet=residuals(mod.mouet, type = "pearson")
N=nrow(mouette)
p=length(coef(mod.mouet))   
sum(res.mouet^2) / (N - p) # forte surdispersion par rapport à un attendu poissonien

# tentative de corriger la surdispersion - modèle négatif binomial
mouet.negbin=glm.nb(comptage~sel+o2+niveau,data=mouette)

par(mfrow=c(2,2))
plot(mouet.negbin) # la déviation à la normalité n'est plus aussi forte qu'avec le modèle poissonnien, et concerne plutôt les 0. Le patron résiduel n'est pas très structuré
plot(qqnorm(residuals(mouet.negbin)),col=col)

res.mouet.nb=residuals(mouet.negbin, type = "pearson")
N=nrow(mouette)
p=length(coef(mouet.negbin))   
sum(res.mouet.nb^2) / (N - p) # on résout bien la surdispersion - on est un peu au dessus de 1, mais peu probable qu'on puisse vraiment mieux faire

# inférence
summary(mouet.negbin) # les comptages de mouettes baissent avec la salinité et l'oxygène dissous, mais ne sont pas influencés par le niveau d'eau

# représentation graphique des effets marginaux et des résidus partiels
p1=plot(ggpredict(mouet.negbin,terms="sel"),residuals=T)+labs(title="Mouette rieuse",xlab="Salinité (g/L)",ylab="Nombre (effet marginal)")
p2=plot(ggpredict(mouet.negbin,terms="o2"),residuals=T)+labs(title="",xlab="O2 (mg/L)",ylab="Nombre (effet marginal)")
p3=plot(ggpredict(mouet.negbin,terms="niveau"),residuals=T)+labs(title="",xlab="Niveau d'eau (cm)",ylab="Nombre (effet marginal)")
plot_grid(p1,p2,p3) # il y a clairement un résidu partiel qui pose problème. Cela ne semble pas affecter le modèle (cf diagnostics), donc nous allons laisser tomber pour l'instant, 
# mais dans un cas réel il faudrait s'en préoccuper : quel est ce comptage? Est-il dans des conditions particulières? 

p1=plot(ggpredict(mouet.negbin,terms="sel"),residuals=T)+labs(title="Mouette rieuse",xlab="Salinité (g/L)",ylab="Nombre (effet marginal)")+ylim(0,1000)
p2=plot(ggpredict(mouet.negbin,terms="o2"),residuals=T)+labs(title="",xlab="O2 (mg/L)",ylab="Nombre (effet marginal)")+ylim(0,1000)
p3=plot(ggpredict(mouet.negbin,terms="niveau"),residuals=T)+labs(title="",xlab="Niveau d'eau (cm)",ylab="Nombre (effet marginal)")+ylim(0,1000)
plot_grid(p1,p2,p3)

# outre les effets des variables, on voit bien que les résidus sont assez dispersés, ce qui suggère une forte variabilité autour de la prédiction
# du modèle. Ce n'est pas grave en soi, mais cela veut dire que nos variables n'expliquent qu'une partie modérée de la variation totale des comptages, qu'on peut quantifier: 

# % déviance expliquée
(951.13-882.30)/951.13 # les variables explicatives représentent 7% de la variation totale des comptages de mouettes. Cela ne veut pas dire qu'elles ont un effet 
# faible ou qu'elles ne sont pas utiles pour expliquer les variations des effectifs de mouettes. Cela veut dire que ces effets sont bien réels, mais qu'ils s'ajoutent
# à de nombreux autres effets qu'il faudrait identifier si on veut espérer prédire avec précision les effectifs de mouettes. Ce n'est pas notre objectif
# ici, on cherche juste à connaître les effets marginaux de ces 3 variables, donc il n'y a rien qui s'oppose à l'inférence et notre conclusion qu'o2 et salinité
# ont un effet négatif sur les mouettes reste valide. On peut imaginer les variables qui expliquent beaucoup la variabilité résiduelle : très probablement la date, 
# éventuellement l'année, mais aussi probablement le bassin sur lequel a lieu chaque comptage. On peut imaginer d'autres facteurs comme les dérangements. 

# incidence ratio rates et intervalles de confiances
irr.mouet.nb=exp(cbind(Estimate = coef(mouet.negbin), confint(mouet.negbin)))
# éventuellement (pas toujours) plus facile à interpréter que les coefficients de pente. On lit que les comptages de mouettes sont multipliés par 0.99 (donc baissent d'1%)
# quand on augmente la salinité d'1g/L, et qu'ils sont multipliés par 0.80 (donc perdent à peu près 1/4) par mg/L d'oxygène dissous.


# est ce que le modèle négatif binomial est mieux ajusté que le modèle poisson (test de Chi² sur les log likelihood)
pchisq(2 * (logLik(mouet.negbin) - logLik(mod.mouet)), df = 1, lower.tail = FALSE) # cette comparaison des log(vraisemblances) nous montre que le modèle negbin est plus
# vraisemblance du modèle poissonnien. Attention, il s'agit ici de vraisemblance statistique (= probabilité que le modèle génère les données). Seule la distribution de la variable
# de réponse a changé entre ces deux modèles, ce qui nous conduit à considérer que la loi négative binomiale est meilleure que la loi de poisson pour modéliser ces comptages surdispersés

# on réajuste le modèle avec variables centrées-réduites pour comparer les effets des paramètres
mouet.negbin.sc=glm.nb(comptage~scale(sel)+scale(o2)+scale(niveau),data=mouette)
summary(mouet.negbin.sc) # cela nous permet de visualiser le plus fort effet de l'oxygène dissout par rapport au sel (2* plus d'effet sur une échelle linéaire) 

# paramètres du modèle de Poisson pour voir l'impact de la surdispersion
irr.mouet=exp(cbind(Estimate = coef(mod.mouet), confint(mod.mouet))) # comparer avec les incidence ratio rates du modèle negbin. En tenant compte de la 
# surdispersion, on élargit les intervalles de confiance, ce qui fait disparaitre un léger effet négatif du niveau d'eau que l'on voyait quand on ne tenait pas 
# compte de la surdispersion --> voyez comment la qualité d'ajustement d'une distribution peut changer vos conclusions 

# on peut tenter une zero inflated Poisson pour tenter de tenir compte encore mieux de la surdispersion des comptages
library(pscl)
zip.mouet=zeroinfl(comptage ~sel+o2+niveau|1,dist="poisson",data=mouette) # dans ce code, on spécifie d'abord les covariables du modèle poissonnien, 
# puis après la | les covariables du modèle binomial.

plot(qqnorm(residuals(zip.mouet))) # clairement, la ZIP ne convient pas. Essayons de mettre des covariables sur le modèle binomial
zip.mouet2=zeroinfl(comptage ~sel+o2+niveau|bassin+dateJ+annee,dist="poisson",data=mouette)
qqnorm(residuals(zip.mouet2)) 
abline(0,1) # pas terrible non plus : le modèle négatif binomial reste le plus approprié

AIC(zip.mouet,mouet.negbin) # le modèle négatif binomial est clairement le plus approprié au regard de la vraisemblance.
