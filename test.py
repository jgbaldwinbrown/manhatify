#!/usr/bin/env python3

import manhatify as mh

def main():
    with open("reflens.bed", "r") as inconn:
        chrlens = mh.get_chrom_lens_from_bed(inconn)
    print(chrlens)
    gdata = mh.get_data_from_bed("answer_dens.bed")
    print(gdata)
    rdata = mh.get_data_from_bed("repanswer_dens.bed")
    print(rdata)
    genedata, chrom_offsets = mh.manhatify(gdata, chrlens, chrom_col = "Scaffold", bp_col = "Start", val_col = "Value", feature = "Genes")
    repdata, chrom_offsets = mh.manhatify(rdata, chrlens, chrom_col = "Scaffold", bp_col = "Start", val_col = "Value", feature = "Repeats")
    combo = mh.combine_data((genedata, repdata))
    print(combo.head())
    print(chrom_offsets)
    mh.plot_manhat(combo, "select_families_dens2.pdf", chrom_offsets, "Value")

if __name__ == "__main__":
    main()
