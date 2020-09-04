
install.packages("circlize")

# load beta coefficient matrix
multivariate_beta<-as.data.frame(fread("data/multivariate_beta.txt"))
multivariate_beta<-multivariate_beta[,-1]
rownames(multivariate_beta)<-colnames(multivariate_beta)
multivariate_beta<-as.matrix(t(multivariate_beta)) # transpose matrx to show linkages by column for circos plotting

library(circlize)
all_states = rownames(multivariate_beta)
n_states = nrow(multivariate_beta)
state_col = c("Age" = "darkgreen",    "Sex" = "darkgreen",
              "BSA" = "darkgreen",  "SBP" = "darkgreen",
              "DBP" = "darkgreen",    "Pulse rate" = "darkgreen",
              "Diabetes" = "darkgreen",     "Smoking" = "darkgreen",
              "Duration of activity" = "darkgreen",     "Medication (n)" = "darkgreen",
              "Assessment centre" = "darkgreen", "HbA1c" = "#377EB8",
              "C-reactive protein" = "#377EB8",  "HDL" = "#377EB8",
              "Glucose" = "#377EB8", "Triglycerides" = "#377EB8","eGFR cystatin" = "#377EB8",
              "Cardiac MRI index" = "#E41A1C",    "PDSRll" = "#E41A1C",
              "PDSRrr" = "#E41A1C",  "Err Global" = "#E41A1C",
              "Ell Global" = "#E41A1C",    "AAo distensibility" = "#E41A1C",
              "DAo distensibility" = "#E41A1C",     "LVSVi" = "#E41A1C",
              "LAVmaxi" = "#E41A1C",     "LAVmini" = "#E41A1C",
              "LVEF" = "#E41A1C","LVCO" = "#E41A1C", 
              "RVSVi" = "#E41A1C")

all_states = names(state_col)

# one for rows and one for columns
state_col2 = c(state_col, state_col)
names(state_col2) = c(rownames(multivariate_beta), colnames(multivariate_beta))

colmat = rep(state_col2[rownames(multivariate_beta)], n_states)
colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)

colmat = paste0(colmat, "A0")
dim(colmat) = dim(multivariate_beta)

circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE) # initialise circos plot

multivariate_chord = chordDiagram(multivariate_beta, col = colmat, grid.col = state_col2,
                       directional = TRUE, annotationTrack = "grid", 
                       big.gap = 10, small.gap = 1) # plot circos for multivariate_beta
circos.clear() # clear this cirsos plot

head(multivariate_chord)
val<-multivariate_chord$value2
p<-which(sign(val)==-1) # make all associations absolute
val[p]<-(-val[p])
multivariate_chord$value1<-val
multivariate_chord$value2<-val
pv<-which(val>=0.4)
pl<-which(multivariate_chord$rn=='PDSRll')
pr<-which(multivariate_chord$rn=='PDSRrr')
pa<-which(multivariate_chord$rn=='LAVmaxi')
pall<-c(pl,pr,pa)
pl<-which(multivariate_chord$cn=='PDSRll')
pr<-which(multivariate_chord$cn=='PDSRrr')
pa<-which(multivariate_chord$cn=='LAVmaxi')
pall2<-c(pl,pr,pa)
pp<-c(pall,pall2)
pval<-c(pp,pv) # include only the positions with beta coefficient > 0.4 apart from the associations 
               # between PDSRll, PDSRrr and LAVmaxi and all other phenotypes.
pval<-unique(pval)
vpall<-val[pval]
multivariate_chord$value1[pval]<-vpall
multivariate_chord$value2[pval]<-vpall

colmat = rep(state_col2[rownames(multivariate_beta)], n_states)
colmat = rgb(t(col2rgb(colmat)), maxColorValue = 255)
dim(colmat) = dim(multivariate_beta)

multivariate_chord$col[-pval]<-paste0(colmat[-pval], "20") # add faint colour in the links with betas < 0.4
head(multivariate_chord)

circos.par(gap.degree=1,canvas.xlim=c(-0.6,0.6), canvas.ylim=c(-1.2,1.2)) # initialise circosplot

chordDiagram(multivariate_chord, col = multivariate_chord$col, grid.col = state_col2,
             directional = TRUE, annotationTrack = c("grid"), link.rank = order(multivariate_chord$col),
             big.gap = 10, small.gap = 1,preAllocateTracks = list(track.height = mm_h(5))) # plot circos for multivariate_chord

## Add sector numbers. The numbers in each sector represent the sum of all coefficients for each specific variable. (Optional)
#
# for(si in get.all.sector.index()) {
#   circos.axis(h = "top", labels.cex = 0.4, sector.index = si, track.index = 2)
# }

# Add names clockwise
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
}, bg.border = NA)

# Add small circular rectangles to represent the proportions of different transitions in each variable
for(i in seq_len(nrow(multivariate_chord))) {
  if(multivariate_chord$value1[i] > 0) {
    circos.rect(multivariate_chord[i, "x1"], -mm_y(1), 
                multivariate_chord[i, "x1"] - abs(multivariate_chord[i, "value1"]), -mm_y(2), 
                col = state_col2[multivariate_chord$cn[i]], border = state_col2[multivariate_chord$cn[i]],
                sector.index = multivariate_chord$rn[i], track.index = 2)
  }
}

circos.clear()

# END