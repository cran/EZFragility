## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.align = "center",
    fig.width = 10,
    fig.height = 8
)

## -----------------------------------------------------------------------------
library(EZFragility)
library(Epoch)
library(ggplot2)
library(ggtext)
library(gsignal)

## -----------------------------------------------------------------------------
dl<-EpochDownloader(progress = FALSE)
names(dl)

## -----------------------------------------------------------------------------
pt01sz1 <- dl$FragilityData_subpt01_1
pt01sz1

pt26sz1 <- dl$FragilityData_subjh103_1
pt26sz1

## -----------------------------------------------------------------------------
butterworthFilter <- function(epoch) {
    order <- 4

    sampling_freq <- metaData(epoch)$samplingRate
    nyquist_freq <- sampling_freq / 2

    lowpass <- 0.5
    highpass <- nyquist_freq * 0.99
    normalized_freqs <- c(lowpass, highpass) / nyquist_freq
    filter_type <- "pass"
    butter_filter <- gsignal::butter(
        n = order, 
        w = normalized_freqs, 
        type = filter_type)

    # Apply filter to epoch data
    mat <- tblData(epoch)
    
    # Apply zero-phase filtering (filtfilt) to each row
    filtered_data <- gsignal::filtfilt(
        filt = butter_filter,
        x = t(mat))

    filtered_data <- t(filtered_data)
    tblData(epoch) <- filtered_data
    
    epoch
}

## -----------------------------------------------------------------------------
pt01sz1Clipped <- pt01sz1 |>
    crop(start=-10, end=10) |>
    butterworthFilter()
pt01sz1Clipped

pt26sz1Clipped <- pt26sz1 |>
    crop(start=-10, end=10) |>
    butterworthFilter()
pt26sz1Clipped

## -----------------------------------------------------------------------------
visualSOZ <- function(epoch, sozNames) {
    p <- plot(epoch)

    elecColor <- rep("black", nrow(epoch))
    elecColor[rownames(epoch) %in% sozNames] <- 'red'
    elecColor <- rev(elecColor) # match the electrode order in the plot

    p + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1)+ 
    theme(axis.text.y = element_markdown(colour = elecColor))
}

## -----------------------------------------------------------------------------
pt01sozName <- rownames(pt01sz1Clipped)[rowData(pt01sz1Clipped)$soz]
pt01Display <- c(pt01sozName, "MLT1", "MLT2", "MLT3", "MLT4")
pt01sz1Reordered <- pt01sz1Clipped[pt01Display, ]
visualSOZ(pt01sz1Reordered, pt01sozName)

## -----------------------------------------------------------------------------
pt26sozName <- rownames(pt26sz1Clipped)[rowData(pt26sz1Clipped)$soz]
excludedElectrodes <- c("RTG29", "RTG30", "RTG31", "RTG32")
pt26sozName <- pt26sozName[!pt26sozName %in% excludedElectrodes]
pt26Display <- c(
    "ABT1", "ABT2",
    pt26sozName[1:16], 
    excludedElectrodes,
    pt26sozName[17:18])
pt26sz1Reordered <- pt26sz1Clipped[pt26Display, ]
visualSOZ(pt26sz1Reordered, pt26sozName)

## -----------------------------------------------------------------------------
library(doSNOW)
# compute fragility

cl <- makeCluster(parallel::detectCores(), type = "SOCK")
registerDoSNOW(cl)

windowNum <- 250
step <- 125
pt01sz1Frag <- calcAdjFrag(epoch = pt01sz1Clipped, window = windowNum, step = step, parallel = TRUE, nSearch=100L, progress = FALSE)

pt26sz1Frag <- calcAdjFrag(epoch = pt26sz1Clipped, window = windowNum, step = step, parallel = TRUE, nSearch=100L, progress = FALSE)

# Stop the parallel backend
stopCluster(cl)

## -----------------------------------------------------------------------------
fragHeatmap <- function(frag, sozNames, ranked=FALSE) {
    startTimes <- frag$startTimes

    indexsz <- which(abs(startTimes)<=0.01)
    elecColor <- rep("black", length(frag$electrodes))
    elecColor[frag$electrodes%in% sozNames] <- 'red'
    elecColor <- rev(elecColor)

    plotFragHeatmap(frag = frag, ranked=ranked) + 
    geom_vline(xintercept = indexsz, color = "black", linetype = "dashed", linewidth = 1) + 
    theme(
        axis.text.y = element_markdown(colour = elecColor)
    )
}

## -----------------------------------------------------------------------------
pt01sz1FragReordered <- pt01sz1Frag[pt01Display]
fragHeatmap(pt01sz1FragReordered, pt01sozName, ranked=FALSE)

## -----------------------------------------------------------------------------
fragHeatmap(pt01sz1FragReordered, pt01sozName, ranked=TRUE)

## -----------------------------------------------------------------------------
pt26sz1FragReordered <- pt26sz1Frag[pt26Display]
fragHeatmap(pt26sz1FragReordered, pt26sozName, ranked=FALSE)

## -----------------------------------------------------------------------------
fragHeatmap(pt26sz1FragReordered, pt26sozName, ranked=TRUE)

## -----------------------------------------------------------------------------
fragDist <- function(frag, sozNames) {
    timeWindows <- frag$startTimes
    timeIdx <- which(timeWindows >= -5 & timeWindows <= 10)
    frag <- frag[, timeIdx]
    plotFragDistribution(frag = frag, groupIndex = sozNames, bandType="SEM", rollingWindow = 1) + 
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1)
}

## -----------------------------------------------------------------------------
fragDist(pt01sz1Frag[pt01Display], pt01sozName)

## -----------------------------------------------------------------------------
fragDist(pt26sz1Frag[pt26Display], pt26sozName)

