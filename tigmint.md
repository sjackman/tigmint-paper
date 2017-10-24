---
title: "Tigmint: Correct Misassemblies Using Linked Reads From Large Molecules"
author: [Shaun D Jackman, Benjamin P Vandervalk, Rene L Warren, Hamid Mohamadi, Justin Chu, Sarah Yeo, Lauren Coombe, Joerg Bohlmann, Steven JM Jones, Inanc Birol]
keywords: [genome, assembly, misassembly correction, scaffolding, 10x, Chromium, linked reads, ABySS, gymnosperm, mitochondrion]
bibliography: tigmint.bib
csl: tigmint.csl
rangeDelim: "&ndash;"
eqnPrefix: "Equation"
figPrefix: "Fig."
tblPrefix: ["Table", "Tables"]
fontsize: 12pt
permalink: /
---

# Abstract

Long-read sequencing technologies have greatly improved assembly contiguity, but at a cost roughly ten times that of short-read sequencing technology. For population studies and when sequencing large genomes, such as conifer genomes and other economically important crop species, this cost may be prohibitive. The 10x Genomics Chromium technology generates linked reads from large DNA molecules at a cost comparable to standard short-read sequencing technologies. Whereas paired-end sequencing gives two reads from a small DNA fragment, linked reads yield roughly a hundred reads from molecules with a typical size of 10 to 100 kilobases. Linked reads indicate which reads were derived from the same DNA molecule, and so should be in close proximity in the underlying genome. Linked reads have been used previously to phase diploid genomes using a reference, *de novo* assemble complex genomes in the gigabase scale, and further scaffold draft assemblies.

In *de novo* sequencing projects, it is challenging yet important to measure the correctness of the resulting assemblies. Linked reads from technologies such as Chromium offer an opportunity to algorithmically address this problem. Here we introduce a software tool, Tigmint, to identify misassemblies using linked reads. The reads are first aligned to the assembly, and the extents of the large DNA molecules are inferred from the alignments of the reads. The physical coverage of the large molecules is more consistent and less prone to coverage dropouts than that of the short read sequencing data. Atypical drops in physical molecule coverage, less than the median minus 1.5 times the inter-quartile range, reveal possible misassemblies. Clipped alignments of the first and last reads of a molecule are used to refine the coordinates of the misassembly with base-pair accuracy.

No software tool currently exists for the specific purpose of identifying misassemblies using linked reads. The tool Long Ranger by 10x Genomics detects structural variants, which is a similar task. It requires however a reference genome assembled in fewer than 500 contigs, whereas a *de novo* assembly is often more fragmented. Tigmint addresses specifically the unaddressed problem of identifying misassemblies using linked reads.

Assemblies of short read sequencing data are easily confounded by repetitive sequence larger than the fragment size of the sequencing library. When the size of a repeat exceeds the library fragment size, the contig comes to an end in the best case, or results in misassembled sequence in the worst case. Tigmint is particularly useful in correcting these misassemblies when the initial assembly of the Illumina paired-end reads did not employ the barcodes of the linked reads, and this rich source of evidence is yet untapped.

Misassemblies not only complicate downstream analyses, but also limit the contiguity of the assembly, when incorrectly assembled sequences prevent joining their adjacent and correctly assembled sequences. To demonstrate the utility of Tigmint, we assemble the six megabase mitochondrial genome of Sitka spruce (*Picea sitchensis*) from 10x Genomics Chromium data using ABySS 2.0, identify and correct misassemblies using Tigmint, and scaffold using ARCS. Tigmint identifies 16 structural misassemblies in this case. After scaffolding with ARCS, the mitochondrial genome is assembled in 12 scaffolds larger than 100 kbp, with an N50 of 493 kbp. We plan to apply this method to assemble the twenty gigabase nuclear genome of Sitka spruce. Chromium reads permit cost-effective assembly of large genomes with high-throughput, short-read sequencing technology, while also providing large-molecule scaffolding data.

The source code of Tigmint is available for download from <https://github.com/bcgsc/tigmint> and is distributed under the GNU GPL v3.0 license.

# Introduction

Genome sequencing reads many short snippets of a genome. Genome assembly attempts to reconstruct the original genome form which these reads were derived. This task is difficult due to gaps and errors in the sequencing data, repetitive sequence in the underlying genome, and heterozygosity, and misassemblies are common. A number of tools exist to correct misassemblies. Pilon [@Walker_2014] maps reads to the assembly and calls variants to corrects small-scale misassemblies. NxRepair [@Murphy_2015] uses Illumina mate-pair sequencing to correct large-scale structural misassemblies.

The linked reads of 10x Genomics Chromium is an Illumina sequencing library prepared from high molecular weight DNA. Short reads derived from the same large molecule are tagged with same 16 nucleotide barcode sequence [@Weisenfeld_2017]. These linked reads have a number of applications. The LongRanger tool of 10x Genomics uses linked reads to map reads to repetitive sequence, phase small variants, and identify structural variants (<https://www.10xgenomics.com/software/>). GROC-SVs [@Spies_2017], NAIBR [@Elyanow_2017], and Topsorter (<https://github.com/hanfang/Topsorter>) identify structural variants using linked reads. ARCS [Yeo_2017], Architect [@Kuleshov_2016], and fragScaff [@Adey_2014] scaffold genome assemblies using linked reads. Supernova [@Weisenfeld_2017] assembles diploid genome sequences using 10x Genomics Chromium data. No tool yet exists for to identify and correct misassemblies using linked reads. We have developed the tool Tigmint for this purpose.

# Methods

## Algorithm

The user provides a draft assembly in FASTA format and the reads in FASTQ format. Tigmint first maps the reads to the reference using BWA-MEM [@Li_2013]. Aligning the reads can alternatively be performed using LongRanger rather than BWA-MEM. The alignments are filtered by mapping quality, alignment score, and number of mismatches to remove poorly aligned reads with the default thresholds $\textrm{MAPQ} > 0$, $\textrm{AS} \geq 100$, and $\textrm{NM} < 5$. Reads with the same barcode that map within 50,000 bp of the adjacent read are grouped into a molecule and assigned a unique numeric molecule identifier. A BED file is constructed, where each record indicates the start and end of one molecule, and records the number of reads that compose that molecule, their median mapping quality, alignment score, and number of mismatches. Molecules shorter than 2000 bp or a user-specified value are removed.

Regions with poor physical molecule coverage indicate potential misassemblies. The depth of molecule coverage at each position is computed from the molecule BED file using Bedtools [@Quinlan_2010]. The median and inter-quartile range (IQR) of the molecule depth of coverage is computed to determine the range of typical physical coverage for the experiment. Regions with coverage less than a threshold specified by the user, and suggested to be the median molecule coverage minus two times the IQR, are flagged as potential misassemblies.

The molecules that span a breakpoint with reads covering that breakpoint will end with a read soft clipped at the position of the breakpoint. These terminal soft clipped reads are used to refine the breakpoint coordinates to base-pair resolution. Four terminal clipped molecules are required by default to identify a breakpoint. The locations of these breakpoints are written to a tab-separated-values (TSV) file along with a summary of the evidence supporting the breakpoint. The sequences of the original draft assembly are split at these breakpoints, producing a FASTA file.

Tigmint will optionally run ARCS [@Yeo_2017] and LINKS [@Warren_2015] at this point to scaffold these corrected sequences into an assembly that is hoped to be more contiguous and correct than the original draft. Tigmint will optionally align the scaffolds to a reference genome, if one is provided, to compute contiguity (NGA50) and correctness (number of breakpoints) metrics of the assembly using ABySS-samtobreak, included with ABySS [@Jackman_2017], before Tigmint, after Tigmint, and after ARCS and LINKS. Each breakpoint identified by ABySS-samtobreak indicates a difference between the assembly and the reference. These breakpoints are composed of both misassemblies and real structural variation between the reference genome and the individual who was sequenced.

## Human data set

We downloaded the ABySS 2.0 [@Jackman_2017] assembly `abyss-2.0/scaffolds.fa` from <http://bit.ly/ncbi-giab-abyss2> of the Genome in a Bottle (GIAB) HG004 Illumina paired-end and mate-pair reads [@Zook_2016]. We downloaded the 10x Genomics Chromium reads for this same individual from <http://bit.ly/giab-hg004-chromium> and used the LongRanger Basic pipeline to extract the barcodes from these reads. We ran Tigmint to correct the ABySS 2.0 assembly of HG004 using these Chromium reads with the command line `tigmint-make depth_threshold=100 draft=abyss2 reads=hg004 ref=GRCh38 G=3088269832`. The median molecule depth of this data is 163, and the IQR is 31, and we set the depth threshold of Tigmint to 100. We provided the reference genome GRCh38 to Tigmint to have it calculate assembly contiguity and correctness metrics. The script to run this analysis is available online at <https://github.com/sjackman/tigmint-data>.

# Results

Breakpoints identified by Tigmint within 1,000 bp of each other were grouped together and counted as a single breakpoint. Tigmint identified 39 breakpoints in the ABySS 2.0 assembly of the GIAB HG004 Illumina paired-end and mate-pair reads using the 10x Genomics Chromium data set from this same individual. The number of breakpoints in this assembly was reduced by 38 breakpoints by using Tigmint, and 97% of the scaffolds modified by Tigmint corrected a breakpoint identified by ABySS-samtobreak. The assembly contiguity (NG50 and NGA50) and correctness (number of breakpoints) metrics before Tigmint, after Tigmint, and after ARCS and LINKS, are shown in @tbl:metrics.

Table: The sequence contiguity and number of breakpoints reported by ABySS-samtobreak when aligned to GRCh38 using BWA-MEM of the ABySS 2.0 assemblies of GIAB HG004. {#tbl:metrics}

| Assembly                   | NG50 (Mbp) | NGA50 (Mbp) | Breakpoints |
| -------------------------- | ---------: | ----------: | ----------: |
| ABySS 2.0                  |       3.49 |        2.97 |       2,717 |
| ABySS 2.0 + Tigmint        |       3.30 |        2.97 |       2,467 |
| ABySS 2.0 + ARCS           |       7.57 |        5.38 |       2,753 |
| ABySS 2.0 + Tigmint + ARCS |      11.54 |        8.98 |       2,493 |

The effect on the precision ($PPV$) and recall ($TPV$) of varying the Depth and Clipped parameters of Tigmint is shown in @fig:precision-recall and @tbl:precision-recall. The assembly is aligned to the reference genome using BWA-MEM. ABySS-samtobreak is used to calculate the number of breakpoints between the assembly and the reference genome. Breakpoints are composed of both misassemblies and true differences, structural variation, between the sequenced individual and the reference genome. Breakpoints due to misassemblies can be corrected, whereas breakpoints due to true structural variation cannot be corrected. The median number of mobile-element insertions, just one class of structural variants, is estimated to be 1,218 per individual [@Sudmant_2015]. For this reason, the sensitivity reported here has an upper bound that is significantly less than perfect.

The number of breakpoints reported by ABySS-samtobreak for the original assembly is the total number of positives, $P = 2,717$, though some represent true structural variation and cannot be corrected. The number of breakpoints identified by Tigmint is the number of predicted positives, $PP = TP + FP$. The number of breakpoints remaining in the assembly after correction by Tigmint is the false negatives, $FN$. The reduction in the number of breakpoints identified by ABySS-samtobreak after the assembly is corrected by Tigmint is the number of true positives, $TP = P - FN$. The number of breakpoints identified by Tigmint minus the true positives is the number of false positives, $FP = PP - TP$. The precision is $PPV = TP / PP$, and the recall is $TPR = TP / P$. The G-score, or Fowlkes–Mallows index, is the geometric mean of precision and recall, $G = \sqrt{PPV \cdot TPR}$. The G-score is maximized when Clipped is 3.

![The effect on precision (PPV) and recall (TPR) of varying the Clipped parameter of Tigmint. The Depth parameter is fixed at 100.](figures/precision-recall.png){#fig:precision-recall}

Table: The effect on precision (PPV), recall (TPR), and G-score of varying the Depth and Clipped parameters of Tigmint. The number of breakpoints detected by ABySS-samtobreak in the uncorrected assembly is $P = 2,717$. The number of breakpoints identified by Tigmint is $PP$. The number of breakpoints remaining in the assembly after correction is $FN$. The reduction in the number of breakpoints is $TP = P - FN$. The number of false positives is $FP = PP - TP$. {#tbl:precision-recall}

| Depth | Clipped | PP   | FN   | TP  | FP   | PPV   | TPR   | G     |
| ----: | ------: | ---: | ---: | --: | ---: | ----: | ----: | ----: |
|   100 |       2 | 1792 | 2467 | 250 | 1542 | 0.140 | 0.092 | 0.113 |
|   100 |       3 |  180 | 2613 | 104 |   76 | 0.578 | 0.038 | 0.149 |
|   100 |       4 |   39 | 2679 |  38 |    1 | 0.974 | 0.014 | 0.117 |
|   100 |       5 |   15 | 2703 |  14 |    1 | 0.933 | 0.005 | 0.069 |

# Discussion

Tigmint uses linked reads to reduce the number of misassemblies in a genome sequence assembly. Scaffolding an assembly that has been so corrected yields an assembly that is both more contiguous and correct than an assembly that has not been corrected.

# References
