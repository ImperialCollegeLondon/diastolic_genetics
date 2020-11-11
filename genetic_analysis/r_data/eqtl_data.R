
to_pull = list(
  
  rs59985551 = data.frame(
    ensembl_id = c("ENSG00000115380"),
    gene = "EFEMP1",
    study = c("GTEx_V8"),
    tissue = c("Thyroid"),
    source = c("api")
  ),
  rs1173727 = data.frame(
    ensembl_id = c("ENSG00000113389"),
    gene = "NPR3",
    study = c("GTEx_V8"),
    tissue = c("Nerve - Tibial"),
    source = c("api")
  ),
  rs35489511 = data.frame(
    ensembl_id = c("ENSG00000157259"),
    gene = "GATAD1",
    study = c("GTEx_V8"),
    tissue = c("Thyroid"),
    source = c("api")
  ),
  rs2275950 = data.frame(
    ensembl_id = c("ENSG00000158006","ENSG00000169442","ENSG00000130695","ENSG00000117632","ENSG00000175087" ,"ENSG00000142675","ENSG00000158014","ENSG00000117676","ENSG00000158022","ENSG00000142669","ENSG00000158008"),
    gene = c("PAFAH2","CD52","CEP85","STMN1","PDIK1L","CNKSR1","SLC30A2","RPS6KA1","TRIM63","SH3BGRL3","EXTL1"),
    study = c("eQTLGen","eQTLGen","eQTLGen","eQTLGen","Fairfax_2014","eQTLGen","GTEx_V8","eQTLGen","TwinsUK","eQTLGen","GENCORD"),
    tissue = c("Blood","Blood","Blood","Blood","monocyte","Blood","Skin - Sun Exposed (Lower leg)","Blood","LCL","Blood","fibroblast"),
    source = c("eqtlgen", "eqtlgen", "eqtlgen", "eqtlgen", "api", "eqtlgen","api","eqtlgen","api","eqtlgen","api")
  ),
  rs11970286 = data.frame(
    ensembl_id = "ENSG00000111860",
    gene = "CEP85L",
    study = c("eQTLGen"),
    tissue = c("Blood"),
    source = c("eqtlgen")
  ),
  rs10261575 = data.frame(
    ensembl_id = "ENSG00000106443",
    gene = "PHF14",
    study = c("eQTLGen"),
    tissue = c("Blood"),
    source = c("eqtlgen")
  ),
  rs11535974 = data.frame(
    ensembl_id = c("ENSG00000259937","ENSG00000139133"),
    gene = c("AC023158.1","ALG10"),
    study = c("HipSci","eQTLGen"),
    tissue = c("iPSC","Blood"),
    source = c("api","eqtlgen")
  ),
  rs499715 = data.frame(
    
  ),
  rs528236848 = data.frame(
    
  ),
  rs9388001 = data.frame(
    ensembl_id = c("ENSG00000111897","ENSG00000025156","ENSG00000152661","ENSG00000146350"),
    gene = c("SERINC1","HSF2","GJA1","TBC1D32"),
    study = c("eQTLGen", "eQTLGen", "GTEx_V8", NA),
    tissue = c("Blood", "Blood", "Nerve - Tibial", NA),
    source = c("eqtlgen", "eqtlgen", "api", "none")
  ),
  rs2234962 = data.frame(
    ensembl_id = c(rep("ENSG00000197771",6), "ENSG00000151929"),
    gene = c(rep("MCMBP",6), "BAG3"),
    study = c("Fairfax_2012","CEDAR","Fairfax_2014","CEDAR","CEDAR","eQTLGen","eQTLGen"),
    tissue = c("B cell","monocyte","monocyte","CD4+ T cell","CD8+ T cell","Blood","Blood"),
    source = c("api", "api", "api", "api", "api", "eqtlgen", "eqtlgen")
  ),
  rs11170519 = data.frame(
    ensembl_id = c("ENSG00000123349","ENSG00000135476","ENSG00000185591","ENSG00000170374","ENSG00000205352","ENSG00000139625","ENSG00000182544"),
    gene = c("PFDN5","ESPL1","SP1","SP7","PRR13","MAP3K12","MFSD5"),
    study = c("GTEx_V8", "GTEx_V8", "eQTLGen", "GTEx_V8", "BLUEPRINT", "eQTLGen", "eQTLGen"),
    tissue = c("Cells - Cultured fibroblasts", "Esophagus - Mucosa", "Blood", "Testis", "monocyte", "Blood", "Blood"),
    source = c("api", "api", "eqtlgen", "api", "api", "eqtlgen", "eqtlgen")
  ),
  rs369533272 = data.frame(
    ensembl_id = c("ENSG00000134779"),
    gene = c("TPGS2"),
    study = c("GEUVADIS"),
    tissue = c("LCL"),
    source = c("api")
  )
  
)

