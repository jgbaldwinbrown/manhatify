#!/usr/bin/env python3

import sys
import pygg
import pandas as pd
import statistics as st
import copy

def get_mid(chromdata_single):
    return(st.mean([chromdata_single["plotpos"].min(), chromdata_single["plotpos"].max()]))

def get_chrom_mids(chromdata, offset):
    chroms = chromdata["Scaffold_number"].unique()
    outframe = pd.DataFrame({"Scaffold_number": chroms})
    mids = []
    for chrom in chroms:
        mids.append(get_mid(chromdata[chromdata["Scaffold_number"].eq(chrom)]))
    outframe["mid"] = mids
    return(outframe)

def get_chrom_lens_from_bed(bedconn):
    chrlens = {}
    for l in bedconn:
        sl = l.rstrip('\n').split('\t')
        chrlens[sl[0]] = int(sl[2])
    return(chrlens)

def get_data_from_bed(bedconn):
    data = pd.read_csv(bedconn, sep="\t", header=None)
    data.columns = ["Scaffold", "Start", "End", "Value"]
    return(data)

def manhatify(indata, chrlens, chrom_col = "chr", bp_col = "bp", val_col = "val", offset = 5e6, feature = "feature"):
    """Take in a data frame and a dictionary of chromosome lengths, then output a tuple with
    two items: an updated data frame with
extra columns containing all of the necessary data for Manhattan plotting, and a data frame
containing chromosomes, rank ordered by size, with their midpoints (in bp) for a manhattan plot"""
    data = copy.deepcopy(indata)
    
    sorted_chroms = sorted(chrlens.items(), key=lambda x: x[1], reverse=True)
    
    chr_lens_dict = {}
    for (chrom, length) in chrlens.items():
        if chrom in data[chrom_col].values:
            chr_lens_dict[chrom] = length
    
    chr_rank_dict = {}
    for rank, (chrom, length) in enumerate(chrlens.items()):
        if chrom in data[chrom_col].values:
            chr_rank_dict[chrom] = rank+1

    cumsum = 0
    cumsum_dict = {}
    for chrom, length in chrlens.items():
        if chrom in data[chrom_col].values:
            cumsum += length
            cumsum_dict[chrom] = cumsum

    data["Scaffold_number"] = data[chrom_col].replace(chr_rank_dict)
    data["chrlen"] = data[chrom_col].replace(chr_lens_dict)
    data["cumsum"] = data[chrom_col].replace(cumsum_dict)
    data["plotpos"] = data["Scaffold_number"] * offset + data["cumsum"] + data[bp_col]
    data["Feature"] = feature
    chrom_mids = get_chrom_mids(data, offset)
    return(data, chrom_mids)

def combine_data(datas):
    return(pd.concat(datas))

def plot_manhat(combo, outname, mids, val_col):
    ggdat = {"data": combo, "outname": outname, "mids": mids}
    
    command = """
    outname = jdata$outname
    mids = jdata$mids
    # labs = c("Gene", "Repeat_content")
    # names(labs) = c("Genes", "Repeats")
    scale = 2.5
    pdf(outname, height=3*4*scale, width=20*scale)
        aplot = ggplot(data = data, aes(plotpos, """ + val_col + """)) +
            geom_point() +
            xlab("Chromosome") + 
            ylab("Features per basepair") + ## y label from qqman::qq
            #scale_color_manual(values = c(gray(0.5), gray(0))) + ## instead of colors, go for gray
            ggtitle("Chromosome-wide feature density") +
            theme_bw() +
            #scale_y_continuous(breaks=seq(from=min(data$""" + val_col + """), to=max(data$""" + val_col + """), length.out=3)) +
            facet_grid(Feature~., scales="free_y") +
            scale_x_continuous(breaks = mids$mid, labels = mids$Scaffold_number) + ## add new x labels 
            theme(text = element_text(size=32))
        print(aplot)
    dev.off()
    
    #   facet_grid(Feature~., scales="free_y", labeller=labeller(Feature=labs)) +
    #   facet_wrap(.~NAME, ncol = 2) +
    #   guides(colour=FALSE) +  ## remove legend
    """
    pygg.ggplot(ggdat, command)
