# load multiplex wormhole functions
source("/Users/maggiehallerud/Desktop/multiplex_wormhole/src/multiplex_wormhole/primerTree_specificity_checks.R")
# set personal API key (for faster processing)
options("ENTREZ_KEY"=ncbi_apikey)


##### STEP 1: RUN PRIMER-TREE FOR EACH PANEL, FILTER RESULTS, & PLOT PRIMER MISMATCHES ######
##### FILTERING DETAILS & JUSTIFICATION IN SECTION 3.1 OF MAIN TEXT & APPENDIX 1 ######
##### IMPORTANT: primerIDs in FASTA headers must be in the format ID.FWD, ID.REV
## black bear
bear <- "BlackBear/Petrou_etal_2025_blackBear_primers.fasta"
FASTA2CSV(bear, gsub(".fasta",".csv", bear))
runPrimerTree(gsub(".fasta",".csv", bear), 
              organisms=c("plants","mammals","birds","fish","bacteria","Homo sapiens"),
              outcsv="BlackBear/Blackbear_primerblast.csv",
              fwd_adapter="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", 
              rev_adapter="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
              THREADS=8,
              PRIMER_SPECIFICITY_DATABASE="nt",
              MAX_AMPLICON_SIZE=300)
bear_pb <- read.csv("BlackBear/Blackbear_primerblast.csv")
unique(bear_pb$family[bear_pb$mismatch_forward<=3 & bear_pb$mismatch_reverse<=3,])
families <- data.frame(family=c("Ursidae","Mustelidae","Felidae","Hominidae",
                                "Muridae","Salicaceae","Poaceae","Canidae",
                                "Suidae","Otariidae","Odobenidae","Phocidae",
                                "Cricetidae",#"Oleaceae",#>3 mismatches
                                "Bovidae","Dasypodidae","Ochotonidae",
                                "Soricidae","Sciuridae","Threskiornithidae",
                                "Leoporidae","Castoridae","Cervidae","Didelphidae",
                                "Gruidae","Cupressaceae"),
                       names=c("Ursidae","Mustelidae","Felidae","humans",
                               "Muridae","Salicaceae","Poaceae","Canidae",
                               "Suidae","Otariidae","Odobenidae","Phocidae",
                               "Cricetidae",#"lilacs",
                               "Bovidae","Dasypodidae","Ochotonidae",
                               "Soricidae","Sciuridae","ibises & spoonbills",
                               "Leoporidae","Castoridae","Cervidae","Didelphidae",
                               "cranes","Cupressaceae"))
bear_pb$name <- families$names[match(bear_pb$family, families$family)]
bear_pb$family <- factor(bear_pb$family, levels=rev(c("Ursidae","Canidae","Felidae",
                                                  "Mustelidae","Odobenidae","Otariidae",
                                                  "Phocidae","Cervidae","Sciuridae","Bovidae",
                                                  "Suidae","Castoridae",
                                                  "Didelphidae","Dasypodidae",
                                                  "Cricetidae","Muridae","Soricidae","Talpidae",
                                                  "Ochotonidae",
                                                  "Poaceae","Oleaceae","Salicaceae",
                                                  "Hominidae")))
plotMismatches(bear_pb[bear_pb$family!="Hominidae",], "family", title="Ursus americanus")+
  theme(axis.text.x=element_blank(), text=element_text(size=10), title=element_text(face="italic"))
ggsave("BlackBear_mismatches.pdf", width=6, height=2.5, dpi=300)


## coyote
coyote <- readxl::read_excel("CALA/ece311266-sup-0001-datas1.xlsx", 1)
write.csv(data.frame(PrimerID=coyote$SNP_id, Sequence=coyote$Primer_sequence_forward),
          "CALA/CALA_primers.csv")
runPrimerTree("CALA/CALA_primers.csv",
              organisms=c("plants","mammals","birds","bacteria","insects"),
              "CALA/CALA_primerblast.csv",
              THREADS=8, fwd_adapter = "",
              rev_adapter = "", PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
cala_pb <- read.csv("CALA/CALA_primerblast.csv")
unique(cala_pb$order)
cala_pb <- cala_pb[cala_pb$order%in%c("Carnivora","Lagomorpha","Rodentia",
                                      "Artiodactyla","Galliformes","Anseriformes",
                                      "Passeriformes","Orthoptera","Coleoptera","Rosales",
                                      "Poales","Ericales","Charadriiformes","Hymenoptera",
                                      "Diptera","Lepidoptera","Perissodactyla"),]
families <- list("Canidae"="dogs",
                 "Felidae"="cats","Mustelidae"="weasels",
                 "Cervidae"="deer","Leporidae"="rabbits, hares",
                 "Suidae"="pigs",
                 "Muridae"="Old world mice and rats",
                 "Balaenopteridae"="rorqual whales",
                 "Balaenidae"="right/bowhead whales",
                 "Physeteridae"="sperm whales",
                 "Zyphiidae"="beaked whales",
                 "Monodontidae"="beluga, narwhal",
                 "Phocoenidae"="porpoises",
                 "Delphinidae"="dolphins",
                 "Eschrichtiidae"="gray whale",
                 "Otariidae"="sea lions, earless seals",
                 "Phocidae"="earless seals",
                 "Odobenidae"="walrus",
                 #"Phasianidae"="turkey, pheasant, quail",
                 #"Anatidae"="waterfowl",
                 "Odontophoridae"="New World quail",
                 "Charadriidae"="plovers",
                 "Scolopacidae"="sandpipers",
                 "Laridae"="gulls",
                 "Corvidae"="crows, jays",
                 "Turdidae"="thrushes",
                 "Icteridae"="blackbirds",
                 "Fringilidae"="finches",
                 "Passerellidae"="sparrows",
                 "Passeridae"="house sparrow",
                 "Hirundinidae"="swallows",
                 "Motacillidae"="pipts",
                 "Tyrannidae"="flycatchers",
                 #"Formicidae"="ants",
                 #"Apidae"="bees",
                 #"Acrididae"="grasshopper",
                 #"Trigonidiidae"="crickets",
                 "Poaceae"="grasses",
                 "Ebenaceae"="persimmons")
cala_pb <- cala_pb[cala_pb$family%in%names(families),]
cala_pb$family <- factor(cala_pb$family, levels=rev(c("Canidae", "Felidae", "Mustelidae", "Cervidae", "Suidae", "Leporidae",
                                             "Muridae", "Phocidae", "Otariidae", "Odobenidae", "Balaenopteridae",
                                             "Balaenidae", "Physeteridae", "Monodontidae", "Phocoenidae", "Delphinidae",
                                             "Eschrichtiidae", "Turdidae", "Icteridae", "Corvidae", "Odontophoridae",
                                             "Passeridae", "Passerellidae", "Charadriidae", "Scolopacidae", "Tyrannidae",
                                             "Motacillidae", "Hirundinidae", "Laridae", "Ebenaceae", "Poaceae")))
plotMismatches(cala_pb, "family", "Canis latrans")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("Coyote_mismatches_300bp.pdf", dpi=300, width=4, height=5)


## great lakes wolves
glw <- "GreatLakesWolves/HERVEY2025_GreatLakeWolves.fasta"
FASTA2CSV(glw, gsub(".fasta",".csv", glw))
runPrimerTree(gsub(".fasta",".csv", glw), 
              organisms=c("Ericaceae","mammals","birds","bacteria","Homo sapiens"),
              "GreatLakesWolves/GLW_primerblast.csv",
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              THREADS=8, PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
glw_pb <- read.csv("GreatLakesWolves/GLW_primerblast.csv")
glw_pb <- glw_pb[glw_pb$family%in%c("Muridae","Canidae","Felidae","Suidae",
                                       "Mustelidae","Ursidae","Otariidae","Talpidae",
                                       "Cricetidae","Balaenopteridae","Bovidae","Phocidae","Sciuridae"),]
glw_pb$family <- factor(glw_pb$family, levels=rev(c("Canidae","Ursidae","Felidae",
                                                      "Mustelidae","Otariidae","Phocidae",
                                                      "Bovidae","Sciuridae","Suidae",
                                                      "Equidae","Balaenopteridae","Cricetidae","Muridae",
                                                      "Talpidae","Phasianidae",
                                                      "Ericaceae","Hominidae")))
plotMismatches(glw_pb[glw_pb$family!="Hominidae",], "family", title="Canis lupus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("GreatLakesWolves_mismatches.pdf", width=7, height=2.5, dpi=300)


## marsh deer
marsh <- "MarshDeer/Wolfenson_etal_2022_southern_marsh_deer.fasta"
FASTA2CSV(marsh, gsub(".fasta",".csv", marsh))
runPrimerTree(gsub(".fasta",".csv", marsh), 
              organisms=c("plants","mammals","bacteria","Homo sapiens"),
              "MarshDeer/MarshDeer_primerblast.csv",
              fwd_adapter="CCCTACACGACGCTCTTCCGATCT", 
              rev_adapter="GTTCAGACGTGTGCTCTTCCGATCT",
              THREADS=8,
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
marsh_pb <- read.csv("MarshDeer/MarshDeer_primerblast.csv")
marsh_pb <- marsh_pb[marsh_pb$family%in%c("Cervidae","Cannabaceae","Hominidae"),]#"Asparagaceae",#>3 mismatches
marsh_pb$family <- factor(marsh_pb$family, levels=rev(c("Cervidae","Cannabaceae","Hominidae")))
plotMismatches(marsh_pb[marsh_pb$family!="Hominidae",], "family", "Blastocerus dichotomus")+
  theme(axis.text.x = element_blank(), title=element_text(face="italic"))
ggsave("MarshDeer_mismatches.pdf", width=8, height=3, dpi=300)##saved as 400x150


## marten
marten <- "Marten/Marten95.fasta"
FASTA2CSV(marten, gsub(".fasta",".csv", marten))
runPrimerTree(gsub(".fasta",".csv", marten), 
              organisms=c("plants","insects","mammals","birds","bacteria","Homo sapiens"),
              "Marten/Marten_primerblast.csv",
              THREADS=8, PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
marten_pb <- read.csv("Marten/Marten_primerblast.csv")
marten_pb <- marten_pb[marten_pb$family%in%c("Hominidae","Mustelidae","Canidae","Felidae","Ursidae",
                                             "Cricetidae","Muridae","Talpidae","Sciuridae","Soricidae",
                                             "Heteromyidae",
                                             "Leporidae","Ochotonidae","Cervidae",
                                             "Phasianidae", "Poaceae"),]
marten_pb$family <- factor(marten_pb$family, levels=rev(c("Hominidae","Mustelidae","Canidae","Felidae","Ursidae",
                                                          "Cricetidae","Muridae","Talpidae","Sciuridae","Soricidae",
                                                          "Heteromyidae",
                                                          "Leporidae","Ochotonidae","Cervidae",
                                                          "Phasianidae", "Poaceae")))
plotMismatches(marten_pb[marten_pb$family!="Hominidae",], "family", "Martes caurina")+
  theme(axis.text.x=element_blank())
ggsave("Marten_mismatches.pdf", width=7, height=2.5, dpi=300)


## NIDGS
nidgs <- "NIDGS/Garrett2024_NIDGS.fasta"
FASTA2CSV(nidgs, gsub(".fasta",".csv", nidgs))
runPrimerTree(gsub(".fasta",".csv", nidgs), 
              organisms=c("plants","insects","mammals","bacteria","Homo sapiens"),
              "NIDGS/NIDGS_primerblast.csv",
              THREADS=8,
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
nidgs_pb <- read.csv("NIDGS/NIDGS_primerblast.csv")
nidgs_pb <- nidgs_pb[nidgs_pb$order%in%c("Rodentia","Poales","Fabales","Rosales",
                                         "Caryophyllales","Brassicales","Diptera",
                                         "Hymenoptera","Orthoptera","Lepidoptera",
                                         "Solanales","Cucurbitales","Asparagales",
                                         "Lamiales","Sapindales","Fagales","Pinales",
                                         "Vitales","Primates"),]
nidgs_pb <- nidgs_pb[-which(nidgs_pb$order=="Primates"&nidgs_pb$family!="Hominidae"),]
nidgs_pb <- nidgs_pb[nidgs_pb$family%in%c("Sciuridae","Poaceae","Fabaceae","Rosaceae",
                                          "Amaranthaceae","Noctuidae","Formicidae","Apidae",
                                          "Trigonidiidae"),]
nidgs_pb$family <- factor(nidgs_pb$family, levels=rev(c("Sciuridae","Poaceae","Fabaceae","Rosaceae",
                                                        "Amaranthaceae","Noctuidae","Formicidae","Apidae",
                                                        "Trigonidiidae")))
plotMismatches(nidgs_pb[nidgs_pb$order!="Primates",], "family", "Urocitellus brunneus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("NIDGS_mismatches.pdf", width=7, height=2, dpi=300)##600x230


## pika
pika <- "Pika/Russelo_pika307.fasta"
FASTA2CSV(pika, gsub(".fasta",".csv", pika))
runPrimerTree(gsub(".fasta",".csv", pika), 
              organisms=c("plants","mammals","bacteria","Homo sapiens"),
              "Pika/Pika_primerblast.csv",
              THREADS=8,
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
pika_pb <- read.csv("Pika/Pika_primerblast.csv")
pika_pb <- pika_pb[pika_pb$product_length<=300,]
pika_pb <- pika_pb[pika_pb$order%in%c("Primates","Lagomorpha",
                                      "Asterales","Poales","Rosales",
                                      "Lamiales","Ericales","Cornales",
                                      "Fabales","Brassicales","Malgiphiales",
                                      "Alismatales","Vitales",
                                      "Sapindales","Selaginellales",
                                      "Fagales","Malvalles","Ranunculales","Myrtales",
                                      "Celastrales"),]
pika_pb <- pika_pb[-which(pika_pb$order=="Primates"&pika_pb$family!="Hominidae"),]
pika_pb$order <- factor(pika_pb$order, levels=rev(c("Primates","Lagomorpha",
                                                      "Asterales","Poales","Rosales",
                                                      "Lamiales","Ericales","Cornales",
                                                      "Fabales","Brassicales","Malgiphiales",
                                                      "Alismatales","Vitales",
                                                      "Sapindales","Selaginellales",
                                                      "Fagales","Malvalles","Ranunculales","Myrtales",
                                                      "Celastrales")))
pika_pb <- pika_pb[pika_pb$family%in%c("Ochotonidae",
                                       "Rosaceae","Poaceae","Salicaceae",
                                       "Brassicaceae","Fabaceae","Cannabaceae"),]
pika_pb$family <- factor(pika_pb$family, levels=rev(c("Ochotonidae","Poaceae","Fabaceae","Brassicaceae",
                                                  "Cannabaceae","Rosaceae","Salicaceae")))
plotMismatches(pika_pb[pika_pb$order!="Primates",], "family", "Ochotona princeps")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("Pika_mismatches2.pdf", dpi=300, width=5, height=2)


## polar bear
polar <- "PolarBear/Haywardetal2022_PolarBear.fasta"
FASTA2CSV(polar, gsub(".fasta",".csv", polar))
runPrimerTree(gsub(".fasta",".csv", polar), 
              organisms=c("plants","mammals","birds","fish","bacteria","Homo sapiens"),
              "PolarBear/Polarbear_primerblast.csv",
              THREADS=8,
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
             PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
polar_pb <- read.csv("PolarBear/Polarbear_primerblast.csv")
families <- list("Ursidae"="bears", 
                 "Canidae"="dogs","Mustelidae"="weasels",
                 #"Felidae"="cats","Procyonidae"="raccoons",#out of range
                 #"Otariidae"="eared seals",#no range overlap
                 "Odobenidae"="walrus","Phocidae"="earless seals",
                 "Monodontidae"="belugas,narwhals",
                 "Balaenopteridae"="rorqual whales",
                 "Balaenidae"="baleen whales",
                 "Eschrichtiidae"="gray whales",
                 "Cervidae"="reindeer,caribou",
                 #"Muridae" = "mice and rats",#undocumented as diet items
                 "Leporidae" = "rabbits and hares",
                 "Cricetidae" = "voles, lemmings, hamsters, New World rats and mice",
                 "Sciuridae" = "squirrels, marmots, ground squirrels",
                 #"Castoridae" = "beavers",#minimal range overlap
                 "Phasianidae" = "pheasants, grouse, ptarmigan, chickens",
                 "Poaceae" = "grasses",
                 "Rosaceae" = "roses, apples, cherries, raspberries, cloudberries",
                 "Fabaceae" = "legumes, peas, beans, clovers",
                 #"Ochotonidae" = "pikas",#minimal habitat overlap
                 "Ericaceae" = "heaths, blueberries, cranberries, crowberries",
                 "Cyperaceae" = "sedges")
polar_pb <- polar_pb[polar_pb$family%in%names(families),]
polar_pb$family <- factor(polar_pb$family, levels=rev(c("Ursidae","Canidae",
                                                        "Felidae","Mustelidae","Procyonidae",
                                                        "Phocidae","Otariidae","Odobenidae",
                                                        "Monodontidae","Balaenopteridae",
                                                        "Balaenidae","Eschrichtiidae","Cervidae",
                                                        "Leporidae","Sciuridae","Castoridae",
                                                        "Cricetidae","Ochotonidae","Muridae",
                                                        #"Equidae",
                                                        "Phasianidae",
                                                        "Ericaceae","Poaceae","Cyperaceae",
                                                        "Rosaceae","Fabaceae"
                                                        #"Hominidae"
                                                        )))
plotMismatches(polar_pb[polar_pb$family!="Hominidae",], "family", "Ursus maritimus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("PolarBear_mismatches.pdf", height=3, width=7, dpi=300)##saved as 700x300


## snow leopard
snow <- "SnowLeopard/Solari_etal_2025_snowLeopard_primers.fasta"
FASTA2CSV(snow, gsub(".fasta",".csv", snow))
runPrimerTree(gsub(".fasta",".csv", snow), 
              organisms=c("plants","mammals","birds","bacteria","Homo sapiens"),
              "SnowLeopard/SnowLeopard_primerblast.csv",
              THREADS=8, fwd_adapter="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", rev_adapter="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
             PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
snow_pb <- read.csv("SnowLeopard/SnowLeopard_primerblast.csv")
families <- list(
  "Hominidae" = "great apes and humans",
  "Felidae" = "cats",
  #"Bovidae" = "cattle, sheep, goats, antelopes",#>3 mismatches
  "Canidae" = "dogs, wolves, foxes",
  #"Poaceae" = "grasses",#>3 mismatches
  "Mustelidae" = "weasels, martens, badgers, otters",
  #"Equidae" = "horses, donkeys, zebras",#>3 mismatches
  "Suidae" = "pigs and wild boar",
  #"Muridae" = "mice and rats",#>3 mismatches
  #"Ursidae" = "bears",#>3 mismatches
  #"Camelidae" = "camels, llamas, alpacas",#>3 mismatches
  #"Phocidae" = "true seals",#>3 mismatches
  "Anatidae" = "ducks, geese, swans")
  #"Corvidae" = "crows, ravens, magpies"#>3 mismatches
snow_pb <- snow_pb[which(snow_pb$family%in%names(families)),]
snow_pb$family <- factor(snow_pb$family, levels=rev(c("Felidae","Mustelidae",
                                                  "Canidae","Ursidae","Phocidae",
                                                  "Bovidae","Camelidae","Equidae",
                                                  "Suidae",
                                                  "Muridae","Castoridae",
                                                  "Corvidae","Anatidae",
                                                  "Poaceae","Hominidae")))
plotMismatches(snow_pb[snow_pb$family!="Hominidae",], "family", "Panthera uncia")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("SnowLeopard_mismatches.pdf", height=1.5, width=5, dpi=300)##save as 450x200


# sitka black-tailed deer (2022 panel)
deer22 <- "SitkaBlackTailedDeer_2022/Russelo_deer244.fasta"
FASTA2CSV(deer22, gsub(".fasta",".csv", deer22))
runPrimerTree(gsub(".fasta",".csv", deer22), 
              organisms=c("plants","mammals","bacteria","Homo sapiens"),
              "SitkaBlackTailedDeer_2022/SitkaDeer22_primerblast.csv",
              THREADS=8, fwd_adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
sitka22 <- read.csv("SitkaBlackTailedDeer_2022/SitkaDeer22_primerblast.csv")
sitka22 <- sitka22[sitka22$order%in%c("Artiodactyla","Fabales","Rosales",
                                      "Lilliales","Caryophyllales","Brassicales",
                                      "Cupressales","Fagales","Primates"),]
sitka22 <- sitka22[-which(sitka22$order=="Primates"&sitka22$family!="Hominidae"),]
sitka22 <- sitka22[-which(sitka22$order=="Artiodactyla"&sitka22$family!="Cervidae"),]
sitka22 <- sitka22[sitka22$family%in%c(#"Fabaceae",
                                       #"Cannabaceae","Fagaceae","Brassicaceae",#>2 mismatches
                                       "Caryophyllaceae",#
                                       "Cupressaceae",
                                       "Cervidae","Hominidae"),]
sitka22$family <- factor(sitka22$family, levels=c("Fabaceae",
                                                  "Cannabaceae","Fagaceae",#>2 mismatches
                                                "Caryophyllaceae","Brassicaceae",
                                                "Cupressaceae",
                                                "Cervidae","Hominidae"))
plotMismatches(sitka22[sitka22$order!="Primates",], "family", "Odocoileus hemionus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("SitkaBTD2022_mismatches.pdf", width=7, height=1.5, dpi=300)


# sitkadeer26
sitkadeer26 <- "SitkaBlackTailedDeer_2026/Stacey2026_deer_reformat.fa"
FASTA2CSV(sitkadeer26, gsub(".fasta",".csv", sitkadeer26))
runPrimerTree(gsub(".fasta",".csv", sitkadeer26), 
              organisms=c("plants","mammals","bacteria","Homo sapiens"),
              "SitkaBlackTailedDeer_2026/SitkaDeer26_primerblast.csv",
              THREADS=8, fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC",
              rev_adapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
sitka26 <- read.csv("SitkaBlackTailedDeer_2026/SitkaDeer26_primerblast.csv")
sitka26 <- sitka26[sitka26$order%in%c("Poales","Malpighiales","Rosales",
                                      "Cucurbitales","Brassicales","Lamiales",
                                      "Solanales","Artiodactyla","Primates"),]
sitka26 <- sitka26[sitka26$family%in%c("Hominidae","Cervidae","Poaceae",
                                       #"Rosaceae","Lamiaceae",#>3 mismatches
                                       "Salicaceae","Brassicaceae"),]
sitka26$family <- factor(sitka26$family, levels=c("Lamiaceae","Brassicaceae","Rosaceae","Salicaceae",
                                                  "Poaceae","Cervidae"))
plotMismatches(sitka26[sitka26$order!="Primates",], "family", "Odocoileus hemionus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("SitkaBTD2026_mismatches.pdf", width=6, height=1.5, dpi=300)


## tiger
tiger <- "Tiger/Natesh_etal_2019_tiger_primers.fasta"
FASTA2CSV(tiger, gsub(".fasta",".csv", tiger))
runPrimerTree(gsub(".fasta",".csv", tiger), 
              organisms=c("plants","mammals","birds","bacteria","Homo sapiens"),
              "Tiger/Tiger_primerblast.csv",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300,
              THREADS=8, fwd_adapter="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", rev_adapter="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")
tiger_pb <- read.csv("Tiger/Tiger_primerblast.csv")
orders <- list(Primates = "primates",
  Carnivora = "carnivores",
  Rodentia = "rodents",
  Perissodactyla = "odd-toed ungulates",
  Artiodactyla = "even-toed ungulates",
  Lagomorpha = "rabbits and hares",
  Dermoptera = "colugos",
  Anseriformes = "waterfowl",
  Pholidota = "pangolins",
  Passeriformes = "perching birds")
tiger_pb <- tiger_pb[tiger_pb$order%in%names(orders),]
families <- list(Hominidae = "great apes",
  Felidae = "cats",
  Canidae = "dogs",
  Ursidae = "bears",
  Mustelidae = "weasels",
  Bovidae = "cattle, antelopes, goats, and sheep",
  Suidae = "pigs",
  Cercopithecidae = "Old World monkeys",
  #Hylobatidae = "gibbons",
  Leporidae = "rabbits and hares",
  Sciuridae = "squirrels",
  Muridae = "mice and rats",
  Anatidae = "ducks, geese, and swans"
  #Manidae = "pangolins",#rare prey
  #Lorisidae = "lorises",#undocumented
  )
tiger_pb <- tiger_pb[tiger_pb$family%in%names(families),]
tiger_pb$family <- factor(tiger_pb$family, levels=rev(c("Felidae","Mustelidae",
                                                        "Canidae","Ursidae",
                                                        "Bovidae","Suidae","Camelidae",
                                                        "Cercopithecidae","Hylobatidae",
                                                        "Leporidae","Sciuridae",
                                                        "Lorisidae","Galagidae","Manidae",
                                                        "Muridae",
                                                        "Anatidae","Hominidae")))
plotMismatches(tiger_pb[tiger_pb$family!="Hominidae",], "family", "Panthera tigris")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("Tiger_mismatches.pdf", height=2, width=6, dpi=300)##saved: 400x200


## wolf microhaplotypes
wolf <- "Wolf_Microhaplotypes/Delomas_etal_2023_wolf_primers.fasta"
FASTA2CSV(wolf, gsub(".fasta",".csv", wolf))
runPrimerTree(gsub(".fasta",".csv", wolf), 
              organisms=c("plants","mammals","birds","bacteria","Homo sapiens"),
              "Wolf_Microhaplotypes/WolfMicrohaps_primerblast.csv",
              THREADS=8, fwd_adapter = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
              rev_adapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300)
microhap_pb <- read.csv("Wolf_Microhaplotypes/WolfMicrohaps_primerblast.csv")
families <- common_names <- list(Hominidae = "great apes",
  Muridae = "mice and rats",
  Felidae = "cats",
  Poaceae = "grasses",
  Canidae = "dogs",
  Anatidae = "ducks, geese, and swans",
  Bovidae = "cattle, antelopes, goats, and sheep",
  Mustelidae = "weasels",
  Columbidae = "pigeons and doves",
  Didelphidae = "opossums",
  #Talpidae = "moles",#>3 mismatches
  #Heteromyidae = "kangaroo rats and pocket mice",
  Castoridae = "beavers",
  Suidae = "pigs",
  Equidae="horses, donkeys",
  Phocidae = "true seals",
  Otariidae = "eared seals",
  Odobenidae="walruses",
  Balaenopteridae="rorquals",
  Balaenidae="right/bowhead whales",
  Phocoenidae="porpoises",
  Ursidae = "bears",
  Cervidae = "deer",
  Cricetidae = "hamsters, voles, and New World rats and mice",
  Rosaceae="roses & blackberries")
  #Fringillidae = "true finches")
microhap_pb <- microhap_pb[microhap_pb$family%in%names(families),]
microhap_pb$family <- factor(microhap_pb$family, rev(c("Canidae","Ursidae","Felidae",
                                                       "Mustelidae","Odobenidae",
                                                       "Phocidae","Otariidae",
                                                       "Cervidae","Castoridae",
                                                       "Suidae","Equidae","Bovidae","Didelphidae",
                                                       "Cricetidae","Muridae","Talpidae",
                                                       "Balaenopteridae","Balaenidae","Phocoenidae",
                                                       "Anatidae","Columbidae",
                                                       "Poaceae","Hominidae")))
plotMismatches(microhap_pb[microhap_pb$family!="Hominidae",], "family", "Canis lupus")+
  theme(axis.text.x=element_blank(), title=element_text(face="italic"))
ggsave("WolfMicrohaps_mismatches.pdf", width=7, height=3, dpi=300)##saved as 700x200



##### COUNT THE NUMBER OF PRIMER PAIRS INVOLVED IN OFF-TARGET AMPLIFICATION #####
countPairs <- function(pb, target=NA, max=3){
  pb$Count <- 0
  pb$Count[which(pb$mismatch_forward<=max & pb$mismatch_reverse<=max)] <- 1
  print(paste("Primer pairs involved in OTA:", 
              length(unique(pb$FWD[which(!(pb$order%in%c("Primates",target)|
                                            pb$family%in%c("Hominidae",target))
                                          & pb$Count==1)]))
  ))
  return(pb)
}#countPairs

### grab count of primer pairs showing OTA w/ <=3 mismatches for each panel
bear_pb <- countPairs(bear_pb, "Ursidae")
cala_pb <- countPairs(cala_pb, "Canidae")
glw_pb <- countPairs(glw_pb, "Canidae")
marsh_pb <- countPairs(marsh_pb, "Artiodactyla")
marten_pb <- countPairs(marten_pb, "Mustelidae")
pika_pb <- countPairs(pika_pb, "Lagomorpha")
nidgs_pb <- countPairs(nidgs_pb, "Rodentia")
polar_pb <- countPairs(polar_pb, "Ursidae")
sitka22 <- countPairs(sitka22, "Artiodactyla")
sitka26 <- countPairs(sitka26, "Artiodactyla")
snow_pb <- countPairs(snow_pb, "Felidae")
tiger_pb <- countPairs(tiger_pb, "Felidae")
microhap_pb <- countPairs(microhap_pb, "Canidae")



##### STEP 2: IDENTIFY PRIMER-BINDING SITES WITHIN THE TARGET SPECIES GENOME #####
runPrimerTree(gsub(".fasta",".csv", bear), 
              organisms=c("Ursus americanus"),
              "BlackBear/BlackBear_primerblast_TARGET.csv",
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              THREADS=8,
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species", 
              MAX_AMPLICON_SIZE=300)
runPrimerTree("CALA/CALA_primers.csv",
              organisms=c("Canis"),#fails with C. latrans
              "CALA/CALA_primerblast_TARGET.csv",
              THREADS=8, fwd_adapter = "", rev_adapter = "",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species", MAX_AMPLICON_SIZE=300)
runPrimerTree(gsub(".fasta",".csv", glw), 
              organisms=c("Canis lupus"),
              "GreatLakesWolves/GLW_primerblast_TARGET.csv",
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              THREADS=8, 
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species", MAX_AMPLICON_SIZE=300)
# marsh deer genome not set up for PRIMER-BLAST
# marten genome not set up for PRIMER-BLAST
# northern Idaho ground squirrel genome not set up for PRIMER-BLAST
runPrimerTree(gsub(".fasta",".csv", pika), 
              organisms=("Ochotona princeps"), 
              "Pika/Pika_primerblast_TARGET.csv",
              THREADS=8,
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species",
             MAX_AMPLICON_SIZE=300)
runPrimerTree(gsub(".fasta",".csv", polar), 
              organisms=c("Ursus maritimus"),
              "PolarBear/Polarbear_primerblast_TARGET.csv",
              THREADS=8,
              fwd_adapter="CTACACGTTCAGAGTTCTACAGTCCGACGATC", 
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", 
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species", 
             MAX_AMPLICON_SIZE=300)
runPrimerTree(gsub(".fasta",".csv", snow), 
              organisms=c("Panthera uncia"),
              "SnowLeopard/SnowLeopard_primerblast_TARGET.csv",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species",
              THREADS=8, fwd_adapter="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", rev_adapter="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
             MAX_AMPLICON_SIZE=300)
runPrimerTree(gsub(".fasta",".csv", deer22), 
              organisms=c("Odocoileus"),#blasts against Odocoileus virginianus- hemionus not set up for PRIMER-BLAST 
              "SitkaBlackTailedDeer_2022/SitkaDeer22_primerblast_TARGET.csv",
              THREADS=8, fwd_adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species")
runPrimerTree(gsub(".fasta",".csv", sitkadeer26), 
              organisms=c("Odocoileus"),#blasts against Odocoileus virginianus- hemionus not set up for PRIMER-BLAST 
              "SitkaBlackTailedDeer_2026/SitkaDeer26_primerblast_TARGET.csv",
              THREADS=8, fwd_adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
              rev_adapter="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species")
runPrimerTree(gsub(".fasta",".csv", tiger), 
              organisms=c("Panthera tigris"),
              "Tiger/Tiger_primerblast_TARGET.csv",
              THREADS=1, fwd_adapter="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", rev_adapter="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
              PRIMER_SPECIFICITY_DATABASE="PRIMERDB/genome_selected_species")
runPrimerTree(gsub(".fasta",".csv", wolf), 
              organisms=c("Canis lupus"),
              "Wolf_Microhaplotypes/WolfMicrohaps_primerblast_TARGET.csv",
              PRIMER_SPECIFICITY_DATABASE="nt", MAX_AMPLICON_SIZE=300,
              THREADS=8, fwd_adapter = "CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
              rev_adapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")



##### STEP 3: CALC GENETIC DISTANCE B/T ON-TARGET & OFF-TARGET SEQS #####
## function to read in genome primer-binding sites, count primer pairs that amplify >2 loci (allows 1 bp mismatch to primers)
readFilterCountGenomeAmps <- function(csv){
  gen <- read.csv(csv)
  # filter
  gen <- gen[which(gen$product_length<=300 & 
                     gen$mismatch_forward<=1 & 
                     gen$mismatch_reverse<=1),]
  # add DNA sequences- reverse complement as needed
  gen$Seq <- Biostrings::DNAStringSet(gen$Sequence)
  gen$Seq[gen$forward_start>gen$reverse_start] <- 
    Biostrings::reverseComplement(gen$Seq[gen$forward_start>gen$reverse_start])
  print(paste("Number of primer pairs with >2 sequences (1 MM allowed): ",sum(table(gen$FWD)>1)))
  return(gen)
}#readFilterCountGenomeAmps

## returns number of primer pairs predicted to amplify >2 sequences for each panel:
uram_gen <- readFilterCountGenomeAmps("BlackBear/Blackbear_primerblast_TARGET.csv")
cala_gen <- readFilterCountGenomeAmps("CALA/CALA_primerblast_TARGET.csv")
glw_gen <- readFilterCountGenomeAmps("GreatLakesWolves/GLW_primerblast_TARGET.csv")
pika_gen <- readFilterCountGenomeAmps("Pika/Pika_primerblast_TARGET.csv")
urma_gen <- readFilterCountGenomeAmps("PolarBear/Polarbear_primerblast_TARGET.csv")
sitka22_gen <- readFilterCountGenomeAmps("SitkaBlackTailedDeer_2022/SitkaDeer22_primerblast_TARGET.csv")
sitka26_gen <- readFilterCountGenomeAmps("SitkaBlackTailedDeer_2026/SitkaDeer26_primerblast_TARGET.csv")
snow_gen <- readFilterCountGenomeAmps("SnowLeopard/SnowLeopard_primerblast_TARGET.csv")
tig_gen <- readFilterCountGenomeAmps("Tiger/Tiger_primerblast_TARGET.csv")
wolf_gen <- readFilterCountGenomeAmps("Wolf_Microhaplotypes/WolfMicrohaps_primerblast_TARGET.csv")


## calcDistToTarget loaded from primerTree_specificity_checks.R
## basic idea: pull all (off-target & on-target) sequences for each primer pair, aligns, 
## then calculates distance of off-target seqs to on-target seq (attempts reverse complement distance if distance >20)
uram_pb <- calcDistToTarget(uram_gen, bear_pb)
# count of primer pairs with offtarget_seqs with <=10 bp difference to the target seq:
length(unique(uram_pb$FWD[which(uram_pb$Dist2Target<=10 & uram_pb$family!="Ursidae")]))
# identify taxonomic families with amplicons that have high similarity to target seqs:
unique(uram_pb$family[which(uram_pb$Dist2Target<=10 & uram_pb$family!="Ursidae")])

cala_pb <- calcDistToTarget(cala_gen, cala_pb)
length(unique(cala_pb$FWD[which(cala_pb$Dist2Target<=10 & cala_pb$family!="Canidae")]))

glw_pb <- calcDistToTarget(glw_gen, glw_pb)
length(unique(glw_pb$FWD[which(glw_pb$Dist2Target<=10 & glw_pb$family!="Canidae")]))

urma_pb <- calcDistToTarget(urma_gen, polar_pb)
length(unique(urma_pb$FWD[which(urma_pb$Dist2Target<=10 & urma_pb$family!="Ursidae")])) 
unique(urma_pb$family[which(urma_pb$Dist2Target<=10 & urma_pb$family!="Ursidae")])

snow_pb2 <- calcDistToTarget(snow_gen, snow_pb)
length(unique(snow_pb2$FWD[which(snow_pb$DistToTarget<=10 & snow_pb$family!="Felidae")]))

tig_pb <- calcDistToTarget(tig_gen, tiger_pb)
unique(tig_pb$FWD[which(tig_pb$Dist2Target<=10 & tig_pb$family!="Felidae")])
unique(tig_pb2$family[which(tig_pb$DistToTarget<=10 & snow_pb$family!="Felidae")])

sitka22_pb <- calcDistToTarget(sitka22_gen, sitka22)
unique(sitka22_pb$FWD[which(sitka22_pb$Dist2Target<=10)])

sitka26_pb <- calcDistToTarget(sitka26_gen, sitka26)
unique(sitka26_pb$FWD[which(sitka26_pb$Dist2Target<=10)])

pika_pb <- calcDistToTarget(pika_gen, pika_pb)
unique(pika_pb$FWD[which(pika_pb$Dist2Target<=10)])

wolf_pb <- calcDistToTarget(wolf_gen, microhap_pb)
length(unique(wolf_pb$FWD[which(wolf_pb$Dist2Target<=10 & wolf_pb$family!="Canidae")]))


### EXAMPLE 1: Primer pair with off-target amplification that's similar to the target amplicon
unique(uram_pb$FWD[which(uram_pb$Dist2Target<=10 & uram_pb$family!="Ursidae")])
rep_region <- uram_pb[uram_pb$FWD=="UrsAm_locus162.1.FWD",]#UrsAm_locus149.1.FWD
rep_gen <- uram_gen[uram_gen$FWD=="UrsAm_locus162.1.FWD",]
rep_merge <- merge(rep_gen, rep_region, all=TRUE)
hist(rep_merge$Dist2Target)
par(mar=c(0,0,1,8))#b, l, t, r
plotAmpliconTrees(rep_merge[rep_merge$mismatch_forward<=1 &rep_merge$mismatch_reverse<=1,],
                  species="Ursus americanus", dG=NA)
write.csv(rep_merge[which(rep_merge$mismatch_forward<=1 & rep_merge$mismatch_reverse<=1 &
                      rep_merge$product_length<=300 & rep_merge$FWD=="UrsAm_locus162.1.FWD"),],
          "BlackBear/UrsAm_locus162_GeneiousExample.csv")
## the above CSV was loaded into Geneious to show a sequence alignment

### EXAMPLE 2: Primer pair that amplifies multiple similar regions within the target species genome
table(pika_gen$FWD)[table(pika_gen$FWD)>1]
pika_gen$SeqNames <- paste(pika_gen$accession, pika_gen$forward_start, sep=":")
plotAmpliconTrees(pika_gen[pika_gen$mismatch_forward<=1&
                             pika_gen$mismatch_reverse<=1&
                             pika_gen$FWD=="Ocp_NRM_3402989_42.1.FWD",],
                  seqnames=pika_gen$SeqNames[pika_gen$mismatch_forward<=1&
                                               pika_gen$mismatch_reverse<=1&
                                               pika_gen$FWD=="Ocp_NRM_3402989_42.1.FWD"],
                  dG=NA, derep=FALSE)#Ocp_NRM_2873624_44.1.FWD, Ocp_NRM_392496_46.1.FWD 


