## ----setup, include = FALSE---------------------------------------------------
library(EZFragility)
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.align = "center",
    fig.width = 6,
    fig.height = 4
)

## ----ictal data---------------------------------------------------------------
data("pt01EcoG")

## create an epoch object
epoch <- Epoch(pt01EcoG)
epoch

## -----------------------------------------------------------------------------
visuIEEGData(epoch  = epoch)

## -----------------------------------------------------------------------------
# The electrode names corresponding to the site of the patient's surgery
sozNames <- attr(pt01EcoG, "sozNames")

## Show the electrodes that are marked as SOZ and additional 4 electrodes
display <- c(sozNames, "MLT1", "MLT2", "MLT3", "MLT4")
visuIEEGData(epoch  = epoch[display])

## Equivalent to: 
## visuIEEGData(epoch  = epoch[display, ])

## constrain to the first 100 time points
visuIEEGData(epoch  = epoch[display, 1:100])

## -----------------------------------------------------------------------------
epochClipped <- truncateTime(epoch, from = -1, to = 0)

visuIEEGData(epoch  = epochClipped)

## -----------------------------------------------------------------------------
## Register a SNOW parallel backend with 4 workers
library(doSNOW)
cl <- makeCluster(4, type = "SOCK")
registerDoSNOW(cl)

windowNum <- 250
step <- 125
pt01Frag <- calcAdjFrag(epoch = epoch, window = windowNum, step = step, parallel = TRUE, nSearch=100L,progress = TRUE)

# Fragility result
pt01Frag

# Stop the parallel backend
stopCluster(cl)

## -----------------------------------------------------------------------------
soz <- estimateSOZ(pt01Frag[ ,pt01Frag$startTimes > 0])
soz

## -----------------------------------------------------------------------------
stats <- fragStat(pt01Frag, soz)
stats

## -----------------------------------------------------------------------------
display <- c(sozNames, "MLT1", "MLT2", "MLT3", "MLT4")
plotFragHeatmap(frag = pt01Frag, sozIndex = sozNames)

## ----out.width="100%"---------------------------------------------------------
plotFragDistribution(frag = pt01Frag, sozIndex = sozNames)

## ----out.width="100%"---------------------------------------------------------
plotFragQuantile(frag = pt01Frag, sozIndex = sozNames)

