---
title: "Assortative mating, sexual selection and their consequences for gene flow in _Littorina_"
author:
  - Samuel Perini:
      institute: tjarno
      correspondence: "yes"
      email: samuel.perini@gu.se
  - Marina Rafajlović:
      institute: tjarno
  - Anja M. Westram:
      institute: ist
  - Kerstin Johannesson:
      institute: tjarno
  - Roger K. Butlin:
      institute: shef
institute:
  - tjarno: Department of Marine Sciences, University of Gothenburg, 40530 Gothenburg, Sweden
  - ist: IST Austria, Am Campus 1, 3400 Klosterneuburg, Austria
  - shef: Department of Animal and Plant Sciences, University of Sheffield, UK, S10 2TN
csl: evolution.xml
bibliography: Lsax_mating.bib
---

The formation of new species requires the evolution of reproductive isolation through the accumulation of barriers to gene flow. Where divergence occurs in allopatry, different barrier effects are automatically associated, but with gene flow these associations need to be created and maintained by selection operating against the effects of recombination [@felsenstein1981; @smadjaandbutlin2011]. One example is the increase in the overall barrier to gene flow resulting from associations between divergent selection and assortative mating [@kirkpatrickandravigne2002; @gavrilets2004; @sachdevaandbarton2017]. If this requires the build-up of linkage disequilibrium among separate sets of loci controlling divergently selected traits, signal traits and preferences, it may be easily opposed by gene flow and recombination [@servedio2009; @smadjaandbutlin2011]. However, some types of traits and forms of assortative mating reduce the number of associations that need to be maintained and so are expected to be more likely to contribute to reproductive isolation. ‘Multiple-effect’ or ‘magic’ traits [@servedio2011; @smadjaandbutlin2011] are traits that contribute to more than one barrier effect. For example, a trait under divergent selection might also function as a mating signal or contribute to mate choice. Assortative mating might depend on a matching rule where signal and preference coincide rather than a preference/trait rule where signal and preference interact [@kopp2018]. In the extreme, there might be only a single trait involved, such as habitat choice or flowering time [“matching rule by a grouping mechanism”; @kopp2018; @servedioandkopp2012]. The ecological trait is then a multiple-effect trait and no other trait is needed to generate assortment. Body size in Gasterosteus sticklebacks [@mckinnonandrundle2002] is an example of a multiple-effect trait where mating is based on phenotypic similarity of a trait under divergent natural selection while wing color-pattern in Heliconius butterflies [@merrill2014; @merrill2019] also being under divergent natural selection, contributes instead to assortative mating primarily through the signal component of a signal-preference system. Assortment can also be driven by the preference component as in the case of cichlids where color sensitivity influences both foraging and mate choice [@seehausen1999].
The evolution of assortative mating, and the barrier to gene flow that it generates, can also be impacted by sexual selection. Assortative mating can occur without variation in mating success among individuals. However, behavioral interactions between males and females that generate assortative mating will often also generate sexual selection. For example, males with intermediate trait values might find mates with common, intermediate preferences more easily than males with extreme values, generating stabilizing sexual selection [@servedio2011; @servedioandhermisson2019]. This stabilizing selection can contribute to reproductive isolation   if the trait optima are different between populations. Sexual selection must be divergent in order to contribute to the ongoing evolution of reproductive isolation but differences in preference between populations may not be enough: if, for example, preferences are less divergent than the traits on which they are based, sexual selection can lead to decreased differentiation between populations after contact [@servedioandboughman2017]. There are still few empirical studies that have demonstrated the extent to which sexual selection contributes to reproductive isolation or its ongoing evolution [@maanandseehausen2011; @servedioandboughman2017].
Whatever the nature of assortative mating and sexual selection, it is important to quantify their contribution to the overall barrier to gene flow during the process of speciation. The contributions of individual barriers can be estimated by breaking down reproductive isolation into its components [@coyneandorr2004 pp. 63-65; @lowry2008; @sobelandchen2014]. In these calculations, the estimate of assortative mating typically comes from comparisons between divergent populations as indices of premating isolation (e.g., Yule’s V and IPSI). In turn, these isolation indices come from experiments where individuals can mate either within their own population or with an individual from a divergent population [e.g., @matsubayashiandkatakura2009]. However, these indices risk over-simplifying the mating pattern and they fail to account for the presence of the intermediate phenotypes that are present whenever reproductive isolation is incomplete [@coyneandorr2004; Irwin 2019].
Hybrid zones provide excellent conditions for quantifying the extent to which gene flow between distinct populations is reduced by divergent natural selection and assortative mating [@hewitt1988]. In contact zones between divergent populations, hybrids can form and display a wide range of trait combinations [@barton1985; @mallet2005]. For example, two locally adapted populations can evolve different trait values for a quantitative trait   but a continuous cline in the trait will typically be maintained across the habitat boundary. Gene exchange will continue but will be impeded, particularly for loci contributing to selected traits and loci closely linked to them. This provides an excellent opportunity to quantify the barrier effects of assortative mating and sexual selection. It has been argued that assortative mating based on clinally-varying traits will generate only a weak barrier to gene flow because individuals that meet one-another in the hybrid zone rarely differ much in trait values, allowing little opportunity for discrimination [@irwin2019]. This logic does not apply to traits with a very simple genetic basis because they are not expected to show a continuous cline across the habitat boundary.  Selection resulting from the reduced fitness of hybrids can, in theory, increase reproductive isolation (reinforcement) but the conditions required are quite stringent [@liou1994; @price2008; but see @servedioandnoor2003]. Both the barrier generated by assortative mating and the likelihood of reinforcement depend on the mechanism of assortment [@kopp2018] and the genetic architecture of the traits involved.
To understand the impact of departures from random mating on the barrier to gene flow in a hybrid zone, it is necessary to quantify the mating pattern. By ‘mating pattern’, we mean the function that predicts the probability of mating, given an encounter between a male and female with specified phenotypes. This might vary across the zone. Given the mating pattern and the distributions of males and female phenotypes, it is possible to predict the strength of assortative mating and sexual selection at any point in the zone. In turn, this can be used to infer the barrier effect in a way that cannot be deduced from interactions between individuals from divergent, parental populations alone. The impacts of assortative mating and sexual selection can also be separated [@servedioandboughman2017] .
Here, we address these issues in the marine snail Littorina saxatilis, combining extensive empirical data from mating experiments with a model-based quantitative description of the mating pattern that we use to infer assortative mating and sexual selection in the field. We then use the mating pattern as an input to computer simulations to study the barrier effects of both assortative mating and sexual selection.
Littorina saxatilis is an intertidal marine snail forming multiple ecotypes, facilitated by low dispersal due to direct development. The Wave and the Crab ecotypes (simply “Wave” and “Crab” in the following) are encountered widely in wave-exposed and crab-rich habitats, respectively, over the species’ North Eastern Atlantic distribution [@johannesson2010; @butlin2014]. Wave individuals live on cliffs, and they have evolved a relatively large foot, thin shell, a bold behavior and small sizes, whereas Crab snails live among boulders, and differ from the Wave snails by a larger, thicker shell with a narrower foot, showing a wary behavior. Trait differences between ecotypes are the result of local adaptation, most likely induced by wave action in the wave-exposed habitat and crab predation in the crab-rich habitat [@johannesson1986; @boulding2017; @lepennec2017]. Many genomic regions potentially involved in the divergence process in L. saxatilis have been identified, including several putative inversions [@westram2018; @faria2019; @morales2019].
Divergent natural selection is a powerful barrier against gene flow between Wave and Crab snail populations but there are also suggestions for other isolating components such as habitat choice and size-assortative mating [@janson1983; @rolan-alvarez1997; @cruz2004; @johannesson2016]. Assortative mating has been investigated in empirical studies both in the field and the laboratory showing that Crab and Wave ecotypes mate assortatively in sympatry [Yule’s V, IPSI and r_i values significantly different from random mating and as high as 0.96; @johannesson1995; @hull1998; @rolan-alvarez1999; @cruz2004; @conde-padin2008] and that female and male sizes in field-collected mating pairs were highly correlated [Pearson correlation coefficients ≥ 0.3; @rolan-alvarez1999; @johannesson1995; @rolan-alvarez2004]. Assortment is accompanied by a component of sexual selection on size that favors large females and small males [@ng2019].   Furthermore, copulation time as well as distances that males followed female trails before mating are longer for similarly sized pairs with the female being on average slightly larger than the male [@hollander2005; @johannesson2008]. Because the average sizes of the ecotypes are very different (adult Crab snails are two to three times larger than adult Wave snails) this generates assortment among ecotypes, with little evidence for effects of traits other than size. Among littorinid snails of various species, males preferentially track and mate females slightly larger than themselves [‘similarity-like’ mechanism plus a constant; @erlandsson1994; @ng2014; @saltin2013; @fernandez-meirama2017; @ng2019] suggesting that this mating pattern is ancestral.
There is strong evidence for the presence of assortative mating by size in L. saxatilis plus the opportunity for sexual selection on size. Thus, size is a multiple-effect trait, under direct divergent selection between the Crab and Wave habitats and also a key trait influencing mating success. However, for the general reasons discussed above, it is unclear to what extent this assortative mating contributes to the barrier to gene flow between the two ecotypes where they meet in natural contact zones. It is also not known whether sexual selection enhances the reproductive barrier in this system. Hence, we asked what the barrier effect of size-assortative mating and sexual selection is in natural contact zones in these snails. First, we quantified the mating probability given encounters between snails with a wide range of sizes and shapes. Second, we used the resulting mating pattern to infer assortative mating and sexual selection across the contact zones between populations of the Crab and Wave ecotypes. Finally, based on the estimates of assortment and sexual selection, we assessed the likely barrier effects of these two components of isolation by performing individual-based computer simulations.

[@ravinet2016] [@rcoreteam2018]  
[@saur1990]. Male mounting position is a reliable proxy for a copulation attempt in L. saxatilis [@hollander2005]. In addition, a positive correlation between mounting duration and the probability that the female received sperm has been observed in other littorinid species [@hollander2018].  
[@lande1981; @gavrilets2004]  
[@glaisher1871]  
[@carpenter2017; @standevteam2018; @rolan2015; @bolkerandteam2017; @derryberry2014]  
[@panova2010; @johannesson2016; @houle1992]  
[@janicke2019; @bartonandbengtsson1986]

# References
