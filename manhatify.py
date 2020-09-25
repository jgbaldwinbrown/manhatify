#!/usr/bin/env python3

import sys
import pygg
import pandas as pd
import statistics as st
import copy

def get_mid(chromdata_single):
    return(st.mean([chromdata_single["plotpos"].min(), chromdata_single["plotpos"].max()]))

def get_chrom_mids(chromdata, offset, chrom_col):
    chroms = chromdata["Scaffold_number"].unique()
    outframe = pd.DataFrame({"Scaffold_number": chroms})
    mids = []
    chrnames = []
    for chrom in chroms:
        mids.append(get_mid(chromdata[chromdata["Scaffold_number"].eq(chrom)]))
        chrnames.append(chromdata[chromdata["Scaffold_number"].eq(chrom)][chrom_col].iloc[0])
    outframe["mid"] = mids
    outframe[chrom_col] = chrnames
    return(outframe)

def get_chrom_lens_from_bed(bedconn):
    chrlens = {}
    for l in bedconn:
        sl = l.rstrip('\n').split('\t')
        chrlens[sl[0]] = int(sl[2])
    return(chrlens)

def get_data_from_bed(bedconn, data_col_name = "Value"):
    data = pd.read_csv(bedconn, sep="\t", header=None)
    data.columns = ["Scaffold", "Start", "End", data_col_name]
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
    if feature:
        data["Feature"] = feature
    chrom_mids = get_chrom_mids(data, offset, chrom_col)
    return(data, chrom_mids)

def combine_data(datas):
    return(pd.concat(datas))

def plot_manhat(combo, outname, mids, val_col, title = "Manhattan plot", dims = (20.0, 10.0), scale = 2, text_size = 32, xname = "Chromosome", yname = "Value", color_col = None, facet_col = None, log = False, named_xticks = False, chrom_col = False):
    ggdat = {
        "data": combo,
        "outname": outname,
        "mids": mids,
        "title": title,
        "dims": dims,
        "scale": scale,
        "text_size": text_size,
        "xname": xname,
        "yname": yname,
        "color_col" : color_col,
        "facet_col" : facet_col
    }
    
    if named_xticks and chrom_col:
        axis_labels = chrom_col
    else:
        axis_labels = "Scaffold_number"
    
    preamble = """
    outname = jdata$outname
    mids = jdata$mids
    scale = jdata$scale
    pdf(outname, height=jdata$dims[2]*scale, width=jdata$dims[1]*scale)
    """
    if color_col:
        aes = "aes(plotpos, " + val_col + ", color = " + color_col + ")"
    else:
        aes = "aes(plotpos, " + val_col + ")"
    
    plot_command = """
        aplot = ggplot(data = data, """ + aes + """) +
            geom_point() +
            xlab(jdata$xname) + 
            ylab(jdata$yname) + ## y label from qqman::qq
            #scale_color_manual(values = c(gray(0.5), gray(0))) + ## instead of colors, go for gray
            ggtitle(jdata$title) +
            theme_bw() +
            scale_x_continuous(breaks = mids$mid, labels = mids$""" + axis_labels + """) + ## add new x labels 
            theme(text = element_text(size=jdata$text_size))"""
    
    if facet_col:
        plot_command = plot_command + """ +
facet_grid(""" + facet_col + """~., scales="free_y")"""
    if log:
        plot_command = plot_command + """ +
scale_y_log10()"""
    plot_command = plot_command + "\n"
    
    ending = """
        print(aplot)
    dev.off()
    
    #   facet_grid(Feature~., scales="free_y", labeller=labeller(Feature=labs)) +
    #   facet_wrap(.~NAME, ncol = 2) +
    #   guides(colour=FALSE) +  ## remove legend
    """
    command = preamble + plot_command + ending
    pygg.ggplot(ggdat, command)
