#---------------------------------------------------------------------------#
### Formation "Analyse de données pour les écologues renforcement" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# TD du mardi : régression multiple
#---------------------------------------------------------------------------#

# packages utiles (NB : aucun n'est indispensable, il ne s'agit que de packages qui facilitent les représentations graphiques)
# certains packages sont re-chargés là où ils sont utilisés afin de faciliter la compréhension du script, mais ce n'est pas une nécessité
library(ggplot2)
library(cowplot)
library(ggpubr)
library(car)
library(jtools)
library(visreg)

# définir le répertoire courant
setwd("D:/certificat_2023/niveau2/supports/donnees")

#---------------------------------------------------------#
### Cas pratique 1 : biométrie chez les manchots adélie ###
#---------------------------------------------------------#

manchots=read.table("manchots_adelie.txt",header=T,sep="\t")

# étude des relations bivariées entre variable de réponse et variables explicatives
p1=ggplot(manchots)+
	aes(x=sexe,y=aile_gauche_moy)+
		geom_boxplot()

p2=ggplot(manchots)+
	aes(x=masse,y=aile_gauche_moy)+
	geom_point()

p3=ggplot(manchots)+
	aes(x=manipulateur,y=aile_gauche_moy)+
	geom_boxplot()

plot_grid(p1,p2,p3)

# étude des colinéarités éventuelles entre variables explicatives : graphiques bivariés
p4=ggplot(manchots)+
	aes(x=sexe,y=masse)+
	geom_boxplot()

p5=ggplot(manchots)+
	aes(x=manipulateur,y=masse)+
	geom_boxplot()

tab.manchot=as.data.frame(table(manchots[,c("manipulateur","sexe")]))

library(ggpubr)
p6=ggballoonplot(tab.manchot,fill="Freq")+
	scale_fill_viridis_c(option = "C")

plot_grid(p4,p5,p6) # il y a un lien masse-sexe 

# on quantifie ce lien au moyen d'une anova - le R² de 0.3 est plutôt rassurant sur le 
# fait que la colinéarité reste faible (il y a un recoupement assez fort), même si la 
# différence de masse entre sexe est statistiquement significative. On fera quand même attention
# à la construction du modèle
summary(lm(masse~sexe,data=manchots))

# distribution de la variable de réponse
hist(manchots$aile_gauche_moy)

#--------------------------#
### modèle de régression ###
#--------------------------#

# modèle
mod.manchots=lm(aile_gauche_moy~masse+sexe+manipulateur,data=manchots)		

# résidus : respect quasi parfait des conditions du modèle linéaire
par(mfrow=c(2,2))
plot(mod.manchots)

# inférence
summary(mod.manchots)

# on détecte des intervalles de confiance assez larges pour la masse et le sexe : ça peut être un effet de la corrélation entre les deux variables. 
# on va creuser en ré-estimant le modèle sans la masse
mod.manchots2=lm(aile_gauche_moy~sexe+manipulateur,data=manchots)		
summary(mod.manchots2)
# puis sans le sexe
mod.manchots3=lm(aile_gauche_moy~masse+manipulateur,data=manchots)		
summary(mod.manchots3)

# puis à sexe constant
mod.manchotsM=lm(aile_gauche_moy~masse+manipulateur,data=manchots,subset=sexe=="M")		
summary(mod.manchotsM)
mod.manchotsF=lm(aile_gauche_moy~masse+manipulateur,data=manchots,subset=sexe=="F")		
summary(mod.manchotsF)

# on note plusieurs points : premièrement, à sexe constant, l'effet de la masse reste sensiblement le même que dans le modèle incluant
# les deux sexes, ce qui suggère que les paramètres ne sont pas biaisés par la colinéarité (le biais est de toute façon une conséquence
# assez rare de la colinéarité, qui a plutôt tendance à augmenter les intervalles de confiance). Deuxièmement, la masse n'a d'effet 
# substantiel dans aucun des deux modèles à sexe constant, donc il n'y a pas vraiment lieu d'attendre un effet fort lorsqu'on inclut
# les deux sexes. 
# dernier point rassurant : le R² du modèle avec les deux variables n'est pas spécialement élevé, suggérant qu'il n'y a pas de surajustement

# si on veut formaliser l'investigation de la colinéarité, on peut se reposer sur les Variance Inflation Factors, qui mesurent l'inflation de la
# variance (= de l'erreur standard) de chaque paramètre du modèle attribuable à la colinéarité entre variables. On considère en général que
# des VIF inférieur à 5 ou 10 (selon les auteurs) impliquent qu'il n'y a pas un gros effet de la colinéarité sur les variances des paramètres.
library(car)
vif(mod.manchots)

# attention néanmoins : un VIF élevé ne signe pas nécessairement un problème, et de manière générale la colinéarité n'est pas nécessairement un 
# problème. La corrélation entre des variables est une réalité biologique qu'il peut être intéressant de reporter dans un modèle. Les situations
# gênantes sont quand (i) deux variables sont colinéaires parce qu'elles représentent la même information, (ii) à l'inverse, lorsque la colinéarité
# n'a aucune base biologique (effets confondants de corrélations fortuites), (iii) quand elles conduisent à un surajustement du modèle.

#------------------------------#
### représentation graphique ###
#------------------------------#

# pour rappel, on représente ici les graphiques de résidus partiels.
# deux solutions, voir d'autres possibilités dans le script de cours

library(ggeffects)
p1=ggpredict(mod.manchots,terms="masse")
p2=ggpredict(mod.manchots,terms="sexe")
p3=ggpredict(mod.manchots,terms="manipulateur")
p1b=plot(p1,residuals=T)
p2b=plot(p2,residuals=T)
p3b=plot(p3,residuals=T)
cowplot::plot_grid(p1b,p2b,p3b)

library(visreg)
par(mfrow=c(2,2))
visreg(mod.manchots,xvar="masse")
visreg(mod.manchots,xvar="sexe")
visreg(mod.manchots,xvar="manipulateur")

# on peut faire un peu de mise en forme
plot(p1,residuals=T)+
	labs(
		x = "Masse (g)", 
		y = "Longueur d'aile (résidus partiels)", 
		title = "")

#--------------------#
### Cas pratique 2 ###
#--------------------#

copro=read.table("coprometrie_orleans.txt",sep="\t",header=T)

# exploration de la variable de réponse (m= masse de fecès de chenilles par échantillon)
par(mfrow=c(2,2))
hist(copro$m) # la distribution des masses est très asymétrique, c'est un grand classique de ce genre d'échantillonnage. On ne pourra pas la modéliser par une loi normale
hist(log(copro$m+1)) # la transformation log n'est pas suffisante (notez l'ajout de la constante 1 pour éviter les log(0) non définis
hist(sqrt(copro$m)) # avec la racine carrée, la distribution est toujours asymétrique, mais c'est nettement moins marqué - à tout le moins, on a bien étalé la variabilité. Cela devrait suffire
copro$sq.mass=sqrt(copro$m) # on prépare d'emblée la variable de réponse transformée à la racine carrée - NB : cette transformation a un coût en termes de facilité d'interprétation

# études des relations bivariées entre la variable de réponse et les variables explicatives
p1=ggplot(copro)+
	aes(x=date_julienne,y=sq.mass)+
		geom_point()

p2=ggplot(copro)+
	aes(x=vol_couronne,y=sq.mass)+
	geom_point()

p3=ggplot(copro)+
	aes(x=diametre_tronc,y=sq.mass)+
	geom_point()

p4=ggplot(copro)+
	aes(x=essence,y=sq.mass)+
	geom_boxplot()

plot_grid(p1,p2,p3,p4) # on voit d'emblée des relations assez fortes se dessiner (avec le volume de la couronne, avec l'essence), et des choses moins claires (tronc, date)

# exploration des relations entre variables explicatives
p5=ggplot(copro)+
	aes(x=vol_couronne,y=diametre_tronc)+
		geom_point()

p6=ggplot(copro)+
	aes(x=essence,y=diametre_tronc)+
	geom_boxplot()

p7=ggplot(copro)+
	aes(x=essence,y=vol_couronne)+
	geom_boxplot()

plot_grid(p5,p6,p7) # on voit de fortes colinéarités se dessiner entre diamètre du tronc et volume de la couronne, et entre ces deux variables et l'essence
# nous n'avons pas regardé les colinéarités entre date julienne et variables liées à l'arbre, parce qu'il n'y a aucune raison biologique pour 
# qu'il y ait une telle corrélation (donc même s'il y en a une, fortuitement, elle ne sera pas dérangeante) et parce que tous les arbres ont été
# échantillonnés à toutes les dates (sauf échec de l'échantillonnage) : ces deux sources de variation sont donc complètement croisées

# on peut pousser plus loin en quantifiant les colinéarités - c'est surtout utile pour la relation diamètre tronc - volume couronne et pour 
# la relation diamètre tronc - essence ; pour volume couronne - essence, la relation est d'emblée si forte qu'on sait que les deux informations seront redondantes
cor.test(copro$vol_couronne,copro$diametre_tronc)# la corrélation est forte, ce qui suggère que les deux informations sont redondante et qu'il sera difficile de les dissocier à l'interprétation du modèle
# il parait judicieux de n'en retenir qu'une - on va donc choisir le diamètre du tronc, qui présente la plus faible relation à l'essence, afin de préserver le plus de variables possible

summary(lm(diametre_tronc~essence,data=copro)) # pour la relation diamètre tronc - essence, on est dans une situation limite où la colinéarité peut ou non poser problème : les deux
# informations sont partiellement redondantes, mais il y a un recoupement suffisamment large pour espérer qu'il reste un peu de variabilité orthogonale. 

# le modèle
mod.copro=lm(sq.mass~date_julienne+essence+diametre_tronc,data=copro)

# vérification des conditions d'application
par(mfrow=c(2,2))
plot(mod.copro) # on voit se dessiner une nette hétéroscédasticité

# le modèle (2)
mod.copro2=lm(sq.mass~date_julienne+essence,data=copro)

par(mfrow=c(2,2))
plot(mod.copro2) #  # ce n'est pas parfait, mais on a bien réduit l'hétéroscédasticité. Il reste aussi une tendance en cloche dans les résidus, qui n'est pas très inquiétante
#  mais suggère que le modèle linéaire que nous avons construit n'est pas la meilleure représentation des données possibles - il y a probablement une contrainte non-linéaire quelque part

# résumé du modèle
summary(mod.copro2) # on voit émerger un fort effet date et un fort effet essence

# que se serait-il passé en laissant les variables colinéaires? 
mod.copro3=lm(m~date_julienne+essence+diametre_tronc,data=copro)
summary(mod.copro3) # l'effet "date" saute sans raison compréhensible

mod.copro4=lm(m~date_julienne+essence+vol_couronne,data=copro)
summary(mod.copro4) # l'effet "date" saute à nouveau, mais aussi l'effet essence

# on peut quantifier l'effet de la colinéarité dans ces différents modèles (Variance Inflation Factor - voir cas pratique 1)
# néanmoins, attention : cette mesure n'est pas infaillible
vif(mod.copro2) # dans ce modèle, il n'y a aucune variable colinéaire - les vif sont faibles
vif(mod.copro3) # dans ce modèle également, malgré un impact évident de la variable "diametre_tronc" qui n'est attribuable à aucune variation interprétable
vif(mod.copro4) # idem ici

mod.copro5=lm(m~date_julienne+essence+vol_couronne+diametre_tronc,data=copro)
vif(mod.copro5) # même dans ce cas extrême, les vif ne mettent en évidence aucun excès de variance, alors qu'on a introduit deux variables très corrélées
summary(mod.copro5)

# message à retenir sur les variance inflation factors et la colinéarité : ce ne sont que des indicateurs qui ne montrent qu'un seul aspect de l'effet de la colinéarité: 
# l'excès de variance sur l'estimation des paramètres, qui n'est ni systématique, ni le seul effet de la colinéarité. Dans le dernier modèle (5), 
# il n'y a ni justification, ni intérêt à introduire vol_couronne et diametre_tronc tant la corrélation entre ces deux variables est forte : il n'y a 
# pas assez de variation orthogonale pour les séparer et par conséquent, ces deux effets sont redondants tant statistiquement que biologiquement, quoiqu'en disent les vif
# vous devez donc faire preuve de jugement critique face à ces indicateurs. 

# représentation graphique des effets

library(ggeffects) # solution avec ggplot
p1=ggpredict(mod.copro2,terms="date_julienne")
p2=ggpredict(mod.copro2,terms="essence")
p1b=plot(p1,residuals=T)
p2b=plot(p2,residuals=T)
cowplot::plot_grid(p1b,p2b) # notez la forte variance autour des deux variables, qui ne les empêche pas d'avoir un effet significatif sur la masse, mais
# qui reflète un bruit élevé malgré un R² pas si mauvais (28%). 

library(visreg) # alternative avec les graphiques R de base
par(mfrow=c(1,2))
visreg(mod.copro2,xvar="date_julienne")
visreg(mod.copro2,xvar="essence")

# Note sur ce cas pratique : ces résultats nécessitent d'être discutés avec précaution du fait du patron quadratique (en cloche) net qui subsiste dans les résidus et
# du reste d'hétéroscédasticité. Les résidus partiels sont intéressants pour essayer de comprendre le problème. Regardez l'agencement des points autour 
# de la droite de l'effet date : ils décrivent une courbure - il semble qu'il y ait un pic de masse quelque part. Nous verrons comment modéliser cela mercredi. 
# Les résidus partiels des essences nous suggèrent que la variabilité des masses en pins est plus faible qu'en chêne : cela pourrait être une source de l'hétéroscédasticité
# que nous n'avons pas réussi à éliminer. Nous verrons comment nous en sortir dans le niveau 3. 
# Un dernier point à discuter : si vous avez bien regardé le jeu de données, vous vous êtes aperçu de deux choses. D'une part, on a échantillonné plusieurs fois chaque
# copromètre - c'est nécessaire pour avoir des mesures comparables à plusieurs dates, mais ça introduit une dépendance entre les individus, violant la condition d'indépendance
# des résidus dans le modèle linéaire, ce qui peut réduire artificiellement les intervalles de confiance. De la même manière, mais moins évident, vous avez peut-être remarqué
# la structure des codes des copromètres qui suggère qu'ils sont regroupés d'une manière ou d'une autre - nous n'avons pas ici les données pour le tester, mais c'est à nouveau une 
# source de non-indépendance des résidus. Nous verrons comment traiter ces effets bloc en niveau 3. 

# On peut donc conclure sur ces résultats, car le modèle linéaire
# est remarquablement robuste et des déviations même assez substantielles à ses conditions d'application sont tolérables - mais si l'interprétation biologique est faisable,
# les conclusions doivent être assorties des avertissements nécessaires.

#----------------------------------------#
#### compléments sur le cas pratique 2 ###
#----------------------------------------#

# Pour ceux qui veulent prendre un peu d'avance et anticiper sur la suite du niveau 2 et le niveau 3, voici quelques solutions. 
# NB : en fin de module 2, on attendra de votre part d'être capable de traiter les deux premiers points de ces compléments et d'identifier
# le problème lié au troisième complément. En fin de module 3, vous serez en mesure de tout traiter. Néanmoins, la détection des problèmes
# est une affaire d'expérience et il est normal de ne pas voir d'emblée certains défauts de modèles, ou de ne pas en percevoir l'importance, 
# quand on commence ses premières analyses. Pratiquez! 

# 1- traitement du patron résiduel. Les graphiques de résidus partiels nous ont suggéré que la date avait un effet quadratique, que nous 
# codons ici. J'ai d'emblée ajouté un terme d'interaction avec l'essence qui s'avère à la fois indispensable et biologiquement intéressant
# (dans la vraie étude, c'est même cette interaction qui nous intéressait spécifiquement, car elle révèle le décalage phénologique entre les 
# chênes et les  pins en termes de disponibilité alimentaire pour les mésanges). En revanche, l'hétéroscédasticité reste assez désastreuse.
library(nlme)
modL1=lm(sq.mass~poly(date_julienne,2)*essence,data=copro)
par(mfrow=c(2,2))
plot(modL1)

# 2- on peut même aller plus loin dans la non-linéarité, ce qui a pour effet de résoudre le problème de l'hétéroscédasticité
library(mgcv)
modL2=gam(sq.mass~s(date_julienne,by=factor(essence)),data=copro)
gam.check(modL2)
plot(predict(modL2),sqrt(residuals(modL2)))
plot(predict(modL2),(residuals(modL2)))
par(mfrow=c(1,2)) # les variations de dates juliennes
plot(modL2)

# 3- afin de tenir compte de la non indépendance des échantillons d'un même copromètre, nous pouvons ajouter un effet aléatoire "copromètre"
# il faudrait aussi tenir compte du regroupement des copromètres dans les différentes parcelles forestières, mais nous n'avons pas les données
# nécessaires ici (dans la réalité nous l'avons bien sûr fait). 
modL3=lme(m~poly(date_julienne,2)*essence,random=list(code_copro=~1),data=copro)
plot(modL3) # l'hétéroscédasticité s'avère particulièrement prononcée
summary(modL3)
modL4=gamm(m~s(date_julienne,by=essence),random=list(code_copro=~1),data=copro)
plot(modL4$lme) # on approche de résidus homoscédastiques - il y a néanmoins peu de données dans les grandes valeurs prédites, mais nous ne pouvons
# rien y changer - les données sont ainsi, on a récolté peu de gros échantillons
par(mfrow=c(1,2))
plot(modL4$gam) # comparez au modèle sans effet aléatoire : la différence n'est pas énorme mais les courbes de prédiction sont aplaties et les intervalles de confiance un peu plus resserrés

# Si on veut pousser encore plus loin, on peut tester l'hypothèse formulée en fin de cas pratique que la variance des masses de fecès dépend de l'essence : 
modL5=gamm(m~s(date_julienne,by=factor(essence)),random=list(code_copro=~1),weights=varIdent(form=~1|essence),data=copro)
plot(modL5$lme) # on rend les résidus encore un peu moins structurés, mais l'amélioration n'est pas flagrante -
par(mfrow=c(1,2))
plot(modL5$gam) # pas de gros changement sur la prédiction de l'effet date
summary(modL5$gam)
summary(modL5$lme) # cette sortie nous confirme qu'il y a plus de variabilité dans les masses côté chêne que côté pins. 
# NB : ce dernier modèle est en limite de compétence attendue en niveau 3. Il commence à être vraiment complexe et implique de savoir ce que vous 
# faites - inutile de l'implémenter si vous n'en comprenez pas chaque terme ou si vous n'avez pas une prédiction biologique claire à chaque variable


