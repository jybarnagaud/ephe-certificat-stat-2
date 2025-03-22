#---------------------------------------------------------------------------#
### Formation "Analyse de données pour les écologues renforcement" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# TD du jeudi : modélisation de structures complexes
#---------------------------------------------------------------------------#
library(ggplot2)
library(cowplot)
library(ggeffects)
library(mgcv)
library(polynom)
library(biostat3)

rm(list=ls())
setwd("C:/Users/jeany/OneDrive/Documents/certificat_2023/niveau2/supports/donnees")

#-----------------------------------------------------------#
### Cas pratique 1 : communautés d'oiseaux dans le Perche ###
#-----------------------------------------------------------#

# charger les données
perche0=read.table("oiseaux_perche.txt",header=T,sep="\t")

# exploration des données
summary(perche0) # attention à la variable essence : il n'y a quasiment que du chêne - il vaut mieux supprimer
# les autres essences, dont les effectifs sont trop faibles pour en tirer une inférence)

perche=droplevels(subset(perche0,essence=="CHE"))
summary(perche)
dim(perche)

# effet massif 
p1=ggplot(perche)+
	aes(x=massif,y=rs)+
		geom_boxplot()

p2=ggplot(perche)+
	aes(x=Hdom,y=rs)+
	geom_point()

p3=ggplot(perche)+
	aes(x=rs)+
	geom_histogram()

plot_grid(p1,p2,p3)

# nb : pour les histogrammes, préférer la fonction graphique de base
hist(perche$rs)

# nb2 : attention : la richesse spécifique suit une distribution de Poisson, même si elle est symétrique
# arguments : c'est une variable discrète bornée sur des valeurs positives et assimilable à un comptage
# la représenter comme une loi normale (non bornée, continue) conduirait à des prédictions fausses

# vérification que moyenne et variance sont proches (définition de la loi de Poisson)
mean(perche$rs)
var(perche$rs)

# d'abord on teste un effet massif
mod.massif=glm(rs~massif,family="poisson",data=perche)
par(mfrow=c(2,2))
plot(mod.massif) # tout va bien

# juste pour la comparaison avec une log transformation de la variable de réponse : regardez le changement des estimations de paramètres
mod.massif.lm=lm(log(rs)~massif,data=perche)

# ensuite on ajoute un effet Hdom au regard de notre hypothèse
mod.mass.hdom=glm(rs~massif+poly(Hdom,2),family="poisson",data=perche)
par(mfrow=c(2,2))
plot(mod.mass.hdom)

# on vérifie la surdispersion mais on se doute qu'il n'y en a pas
res.mass.hdom=residuals(mod.mass.hdom, type = "pearson")
N=nrow(perche)
p=length(coef(mod.mass.hdom))   
sum(res.mass.hdom^2) / (N - p) # les données sont même sous dispersées

# résultats de ce modèle : 

p1=ggemmeans(mod.mass.hdom,terms="massif") # avec ggemmeans, les effets massif sont estimés à Hdom moyen
plot(p1,residuals=T)

p1b=ggpredict(mod.mass.hdom,terms="massif") # avec ggpredict, les effets massif sont estimés à Hdom=0
plot(p1b,residuals=T)
# ici, la différence entre ggemmeans et ggpredict n'est pas très grande, mais elle peut avoir un gros effet sur des modèles
# à plus de variables ou quand les effets sont plus forts. Il ne s'agit que d'une modification de l'apparence des résultats, 
# pas des résultats eux-mêmes : vous pouvez donc choisir la représentation la plus pertinente par rapport à votre objectif. 
# par exemple, ici, je considère qu'un effet massif à Hdom=0 n'a pas beaucoup de sens vu que les parcelles coupées ne sont 
# pas notre principal intérêt

p2=ggemmeans(mod.mass.hdom,terms="Hdom") # effet hauteur dominante pour un massif moyen (avec ggpredict, ce serait l'effet hauteur pour le massif de référence)
plot(p2,residuals=T)

# incidence rate ratio
eform(mod.mass.hdom) # les incidence rate ratios pour des termes polynomiaux sont particulièrement ininterprétables --> on est obligés de passer par une sortie graphique

# est-ce que le modèle à polynome est plus pertinent que le modèle à effet linéaire?
mod.mass.hdom.lin=glm(rs~massif+Hdom,family="poisson",data=perche)
AIC(mod.mass.hdom.lin,mod.mass.hdom) # oui

# NB cette étape n'est pas nécessaire. Si votre hypothèse implique un terme quadratique, alors le terme quadratique suffit. S'il est estimé 
# proche de 0, cela reviendra à un effet linéaire. La comparaison par AIC n'est utile que si une de vos question implique de comparer deux 
# hypothèses concurrentes, une hypothèse linéaire à une hypothèse quadratique, chacune fondée sur une argumentation biologique. Elle peut 
# être aussi utile si vous avez absolument besoin de parcimonie (par ex. visée prédictive mais avec relativement peu de données). Nous montrons
# cette comparaison ici à titre d'illustration uniquement. 

# deuxième hypothèse : on teste l'interaction massif * effet quadratique
mod.inter=glm(rs~massif*poly(Hdom,2),family="poisson",data=perche)
par(mfrow=c(2,2))
plot(mod.inter)
summary(mod.inter)

p3=ggemmeans(mod.inter,terms=c("Hdom","massif"))
plot(p3,residuals=T)
plot(p3,residuals=T,facet=T) # plus lisible en séparant les massifs

# effet interaction ou pas?
AIC(mod.inter,mod.mass.hdom) # l'interaction apporte de l'information au modèle = les données supporte l'idée que les relations
# richesse - hauteur du peuplement diffèrent d'un massif à l'autre. A vous de voir si c'est pertinent biologiquement. 

# on aurait aussi pu traiter la question en gam
mod.mass.hdom.gam=gam(rs~s(Hdom,by=factor(massif)),family="poisson",data=perche)
gam.check(mod.mass.hdom.gam)
par(mfrow=c(2,2))
plot(mod.mass.hdom.gam)

# quelle est la hauteur de peuplement pour laquelle la richesse est maximale / minimale
# vu qu'il y a un terme d'interaction, on va devoir trouver une solution par massif
y=function(x){coef(mod.inter)[1]+
									coef(mod.inter)[5]*x+
												coef(mod.inter)[6]*x*x} # massif de Bellême
bell.h=subset(perche,massif=="Belleme")$Hdom
xmin.bell=optimize(y, interval=range(bell.h), maximum=F) # attention : optimiser sur la range de hauteurs du massif concerné
xmin.bell # hauteur d'arbres pour la richesse minimum
xmax.bell=optimize(y, interval=range(bell.h), maximum=T)
xmax.bell # hauteur d'arbres pour la richesse maximum (!! min et max ne s'entendent que sur la gamme de hauteurs échantillonnée, 
# sans préjuger de ce qui se passerait pour des arbres plus hauts / bas)

# même chose pour Bourse
bour.h=subset(perche,massif=="Bourse")$Hdom
y=function(x){coef(mod.inter)[1]+coef(mod.inter)[2]+
		(coef(mod.inter)[5]+coef(mod.inter)[7])*x+
			(coef(mod.inter)[6]+coef(mod.inter)[10])*x*x}

xmin.bour=optimize(y, interval=range(bour.h), maximum=F)
xmin.bour # hauteur d'arbres pour la richesse minimum
xmax.bour=optimize(y, interval=range(bour.h), maximum=T)
xmax.bour 

# même chose pour Moulins
moul.h=subset(perche,massif=="Moulins")$Hdom
y=function(x){coef(mod.inter)[1]+coef(mod.inter)[3]+
		(coef(mod.inter)[5]+coef(mod.inter)[8])*x+
		(coef(mod.inter)[6]+coef(mod.inter)[11])*x*x}

xmin.moul=optimize(y, interval=range(moul.h), maximum=F)
xmin.moul # hauteur d'arbres pour la richesse minimum
xmax.moul=optimize(y, interval=range(moul.h), maximum=T)
xmax.moul 

# même chose pour Reno
rv.h=subset(perche,massif=="Reno_Valdieu")$Hdom
y=function(x){coef(mod.inter)[1]+coef(mod.inter)[4]+
		(coef(mod.inter)[5]+coef(mod.inter)[9])*x+
		(coef(mod.inter)[6]+coef(mod.inter)[12])*x*x}

xmin.rv=optimize(y, interval=range(rv.h), maximum=F)
xmin.rv # hauteur d'arbres pour la richesse minimum
xmax.rv=optimize(y, interval=range(rv.h), maximum=T)
xmax.rv 

minmax=data.frame(foret=c("Bellême","Bourse","Moulins","Reno Valdieu"),
mini=c(xmin.bell$minimum,xmin.bour$minimum,xmin.moul$minimum,xmin.rv$minimum),
	maxi=c(xmax.bell$maximum,xmax.bour$maximum,xmax.moul$maximum,xmax.rv$maximum))

x11()
ggplot(minmax)+
		geom_segment(aes(x=foret,xend=foret,y=mini,yend=maxi))
# la gamme de variation de hauteur prédite pour la richesse spécifique est à peu près similaire d'un massif à l'autre - un peu plus réduite à Moulins
# cela suggère que malgré le fait que l'interaction massif*hauteur soit statistiquement significative, elle n'a pas un effet majeur sur
# le gradient de richesse spécifique. De fait, dans le dernier graphique de résultats, on voyait déjà se dessiner des courbes assez similaire
# Ce résultat illustre bien la différence entre significativité statistique et biologique d'un résultat : l'interaction est certes différente 
# de 0, mais sa présence dans le modèle ne change rien en termes d'interprétation. 

#---------------------------------------------------------#
### Cas pratique 2 : phénologie horaire des chiroptères ###
#---------------------------------------------------------#

# prepdata (à garder ailleurs) 
chiro=read.table("activite_chiropteres.txt",header=T,sep="\t")

## quelques manipulations préliminaires pour faciliter les analyses (NB : le jeu de données dont vous disposez a déjà intégré ces manipulations)

# passer la date ("Date_nuit") en date julienne, ces commandes sont donc inactivées
	#library(lubridate)
	#chiro$julian_date=yday(lubridate::dmy(chiro$Date_nuit))

# passer les colonnes heures-minutes en heure décimale. De plus, on choisit de centrer la journée autour de minuit pour éviter de scinder
# la nuit en deux groupes horaires
	#chiro$heure.dec=round(chiro$Heure+(chiro$Minute/60),2)
	#chiro$heure.rel=chiro$heure.dec
	#chiro[which(chiro$heure.dec>12 & chiro$heure.dec<=24 ),"heure.rel"]=chiro[which(chiro$heure.dec>12 & chiro$heure.dec<=24 ),"heure.rel"]-24

## filtrage des données

# L'enjeu de ce cas pratique est de sélectionner les bonnes données. Il faut être représentatif, tout en évitant qu'une surcharge
# de données sur une certaine espèce, période horaire ou saison ne vienne polluer l'analyse. 

# répartition des données par période
hist(chiro$julian_date) 
summary(factor(chiro$Mois)) # les mois de juin et juillet génèrent clairement le plus grand nombre de données : on ne garde que ces deux mois
chiro2=droplevels(subset(chiro,Mois%in%c(6,7))) # penser à mettre à jour les variables facteurs avec droplevels. On note au passage que quelques espèces
# vont disparaître à cause de ce choix phénologique. Si ce n'est pas acceptable pour vous, privilégiez le critère "espèce" par rapport au critère "phéno"

# répartition des espèces
summary(chiro2$Espece_validee) # l'hétérogénéité interspécifique (en nb de données) est très forte. Il y a des espèces sur-représentées, comme
# la pipistrelle, et des espèces sous représentées, comme l'oreillard roux (Pleaur). Il va falloir faire un choix. Si on retient les deux pipistrelles
# communes, elles vont prendre toute la variance de notre modèle, créer une très forte surdispersion qu'il nous sera difficile de traiter, et écraser
# les résultats (= le modèle sera complètement dirigé par ces deux espèces). Il vaudrait mieux les modéliser à part. Inversement, vu la complexité du 
# modèle, on va préférer conserver assez de données pour que les estimations soient robustes - on met un seuil arbitraire à 100 données par espèce, 
# ce qui nous conduit à retenir 5 espèces.

chiro3=droplevels(subset(chiro2,Espece_validee%in%c("Barbar","Myosp","Nyclei","Plesp","Rhihip")))

## exploration des données

ggplot(chiro3)+
	aes(x=heure.rel,y=NbCris)+
	geom_point()+
		facet_wrap(~Espece_validee,ncol=2)

ggplot(chiro3)+
	aes(x=julian_date,y=NbCris)+
	geom_point()+
	facet_wrap(~Espece_validee,ncol=2) # attention, 2 campagnes distinctes --> passer julian_date en facteur

chiro3$Mois=factor(chiro3$Mois)

ggplot(chiro3)+
	aes(x=Mois,y=NbCris)+
	geom_boxplot()+
	facet_wrap(~Espece_validee,ncol=2)

hist(chiro3$NbCris)
mean(chiro3$NbCris)
var(chiro3$NbCris) # le nombre de cris est surdispersé par rapport à une loi de Poisson

# on essaie quand même un modèle de Poisson, mais il sera probablement inadéquat
chiro.glm=glm(NbCris~(Mois+poly(heure.rel,2))*Espece_validee,data=chiro3,family=poisson)
# on contrôle d'emblée la surdispersion
res.glm=residuals(chiro.glm, type = "pearson")
N=nrow(chiro3)
k=length(coef(chiro.glm))   
sum(res.glm^2) / (N - k) # les résidus sont surdispersés, mais pas excessivement -> on tente un glm négatif-binomial

# GLM negbin
library(arm)
chiro.glm1=glm.nb(NbCris~(Mois+poly(heure.rel,2))*Espece_validee,data=chiro3)

# surdispersion?
res.glm1=residuals(chiro.glm1, type = "pearson")
N=nrow(chiro3)
k=length(coef(chiro.glm1))   
sum(res.glm1^2) / (N - k) # on peut en rester à un GLM négatif binomial, la surdispersion résiduelle a presque disparu

par(mfrow=c(2,2))
plot(chiro.glm1) # il y a un peu d'hétéroscédasticité, probablement lié au déficit de points dans les faibles valeurs prédites.
# assez classiquement, ce genre de patron est dû à un groupe de données bien identifiable (par ex une espèce ou un habitat rares). 
# ici, il semble logique d'explorer la possibilité qu'une espèce soit responsable de ce patron.

res=data.frame(esp=chiro3$Espece_validee,resid=res.glm,pred=predict(chiro.glm1))

ggplot(res)+
	aes(x=pred,y=resid,colour=esp)+
	geom_point() # sans grande surprise la queue de distribution est due à une seule espèce, Nyclei (Noctule de Leisler). On peut l'enlever et refaire le modèle
# NOTE : dans votre M&M, vous indiquerez clairement que vous avez supprimé cette espèce à cause de la surdispersion malgré le fait qu'elle
# rentre dans le critère des 100 données. Vous indiquerez aussi que vous avez tenté l'analyse avec vs sans l'espèce, et que cela ne change
# pas grand chose aux résultats (voir plus bas). Vous pouvez mettre l'analyse avec Nyclei en annexe de vos résultats, ou au minimum la tenir à la 
# disposition des relecteurs s'ils la demandent

# modèle sans Nyclei
chiro4=droplevels(subset(chiro3,Espece_validee!="Nyclei"))
chiro.glm2=glm.nb(NbCris~(Mois+poly(heure.rel,2))*Espece_validee,data=chiro4)

par(mfrow=c(2,2))
plot(chiro.glm2) # on a fait disparaitre la surdispersion 

#est-ce que ça change quelque chose aux estimations pour les espèces restantes? on compare les coefficients de régression du modèle avec/sans Nyclei
coef.1=coef(chiro.glm1)
coef.2=coef(chiro.glm2)[names(coef.1)]
plot(coef.1,coef.2,type="n")
text(coef.1,coef.2,labels=names(coef.1),cex=0.5)
abline(0,1) # les coefficients sont bien corrélés même si quelques uns diffèrent un peu
# (NB ce n'est pas étonnant qu'il y ait peu de différence, car les espèces sont indépendantes dans le modèle)
# afin de produire une analyse parfaitement propre, on va tout de même poursuivre sans Nyclei 

# on va tester la nécessité des termes d'interaction : cette sélection vient du fait que les interactions représentent notre
# questionnement sur les différences interspécifiques dans les phénologies mensuelles et horaires
chiro.glm3=glm.nb(NbCris~(poly(julian_date,2)+poly(heure.rel,2))+Espece_validee,data=chiro4)
chiro.glm4=glm.nb(NbCris~poly(julian_date,2)*Espece_validee+poly(heure.rel,2),data=chiro4)
chiro.glm5=glm.nb(NbCris~poly(julian_date,2)+Espece_validee*poly(heure.rel,2),data=chiro4)
AIC(chiro.glm2,chiro.glm3,chiro.glm4,chiro.glm5) # le meilleur modèle maintient les deux interactions, mais on est près du seuil de 2 unités d'AIC
# vu nos objectifs, le fait de garder ou non l'interaction n'aura pas une très grande importance, mais si la sélection de modèle était 
# le coeur de l'analyse, on serait plus prudents à ce sujet

# résultats
summary(chiro.glm2) # ininterprétable, trop de paramètres. On va interpréter graphiquement

p1=ggemmeans(chiro.glm2,terms=c("heure.rel","Espece_validee"))
plot(p1,residuals=T,facet=T) # il faut identifier le résidu extrême pour Barbar, car il peut avoir une incidence sur les paramètres de régression

# recherche du résidu extrême
index.barbar=which(chiro4$Espece_validee=="Barbar")
resid.barbar=residuals(chiro.glm2)[index.barbar]
chiro4[names(which.max(resid.barbar)),] # cette ligne correspond au résidu aberrant
barb=subset(chiro4,Espece_validee=="Barbar")
hist(barb$NbCris) # on voit ce point de données seul à droite

par(mfrow=c(1,2)) # on va visualiser où se trouve ce point extrême dans les deux gradients qui nous concernent
plot(chiro4$julian_date,chiro4$NbCris) 
points(chiro4[names(which.max(resid.barbar)),"julian_date"],chiro4[names(which.max(resid.barbar)),"NbCris"],col="red",pch=21,bg="red")

plot(chiro4$heure.rel,chiro4$NbCris)
points(chiro4[names(which.max(resid.barbar)),"heure.rel"],chiro4[names(which.max(resid.barbar)),"NbCris"],col="red",pch=21,bg="red")

# on va refaire le modèle sans la valeur très élevée, qui est en plein dans le gradient des 2 variables étudiées mais
# systématiquement à part du nuage de points : cela peut avoir de nombreuses origines = erreur d'identification, problème d'enregistreur,
# individu tournant autour de l'enregistreur, afflux d'individus pour une raison inconnue - en tout cas, ce n'est pas lié à l'une des 
# deux variables qui nous intéressent. 
# là encore, attention à la rédaction des résultats. Vous avez deux possibilités : (1) dire clairement que vous retirez un point qui 
# vous semble aberrant du point de vue de la variable de réponse, mais n'est pas extrême sur les variables explicatives. Dans ce cas, 
# vous devez faire l'analyse avec/sans le point et présenter l'une en tenant l'autre à disposition du lecteur (annexe, ou disponible sur 
# demande). (2) choisir de conserver le point, dire que vous avez tenté l'analyse sans (en tenant à nouveau les résultats à disposition), 
# mais en étant clair sur le fait qu'il s'agit d'un résidu anormal. 
# Dans tous les cas, nous ne nous autorisons au retrait de ce point que parce que (1) nous avons de bonnes raisons de croire qu'il s'agit
# d'une forme d'erreur de mesure, (2) il n'est pas extrême au regard des covariables et donc porte peu d'information, (3) nous sommes 
# transparents là dessus auprès du lecteur, (4) son retrait ne change pratiquement pas les résultats. Si ces 4 conditions ne sont pas toutes
# simultanément remplies, vous ne pouvez pas vous permettre de supprimer un point de données. 

# analyse du modèle
chiro5=chiro4[-which(rownames(chiro4)==names(which.max(resid.barbar))),]
chiro.glm2b=glm.nb(NbCris~(Mois+poly(heure.rel,2))*Espece_validee,data=chiro5)
chiro.glm3b=glm.nb(NbCris~(poly(julian_date,2)+poly(heure.rel,2))+Espece_validee,data=chiro5)
chiro.glm4b=glm.nb(NbCris~poly(julian_date,2)*Espece_validee+poly(heure.rel,2),data=chiro5)
chiro.glm5b=glm.nb(NbCris~poly(julian_date,2)+Espece_validee*poly(heure.rel,2),data=chiro5)
AIC(chiro.glm2b,chiro.glm3b,chiro.glm4b,chiro.glm5b) # les résultats ne changent pas

summary(chiro.glm2b)

coef.1=coef(chiro.glm2)
coef.2=coef(chiro.glm2b)[names(coef.1)]
plot(coef.1,coef.2,type="n")
text(coef.1,coef.2,labels=names(coef.1),cex=0.5)
abline(0,1) # la valeurs aberrante ne change absolument rien aux estimations, même pour Barbar

## résultats et inférence sur le modèle retenu en dernière instance
p2=ggemmeans(chiro.glm2b,terms=c("heure.rel","Espece_validee"))
x11()
plot(p2,residuals=T,facet=T) # on n'y voit pas grand chose à cause de l'échelle de Barbar, beaucoup plus étalée que les autres + forte variabilité 
plot(p2,residuals=F,facet=T) # sans les résidus (mais avoir conscience que les données s'étalent beaucoup vers le haut : si vous choisissez
# de mettre ce graphique dans vos résultats de publication, vous devez mettre aussi l'autre en annexe pour bien montrer la variabilité des données

p3=ggemmeans(chiro.glm2b,terms=c("Mois","Espece_validee"))
plot(p3,residuals=T,facet=T)  
plot(p3,residuals=F,facet=T) 

# Il y a relativement peu de variation d'activité entre juin et juillet pour les 4 espèces considérées, même si Myosp tend à être un 
# peu plus actif en juillet. Myosp est le plus actif (ou du moins, le plus détecté) - penser cependant aux données brutes : l'étalement 
# des données de Barbar est plus fort que pour les autres espèces. 

# les phénologies horaires nous montrent, là encore avec beaucoup de variation résiduelle, que Barbar et dans une moindre mesure Rhihip
# sont plus actives en début et fin de nuit, alors que Myosp et Plesp sont plus actif vers 1h30 du matin.

# Vous devrez comparer ces résultats à la littérature afin d'en vérifier la cohérence biologique. Un défaut de cohérence peut être dû à 
# une situation exceptionnelle sur le lieu de votre échantillonnage, un biais de mesure, une source d'hétérogénéité cachée, ou un problème
# quelconque dans les données. 

