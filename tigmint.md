---
title: "Tigmint: Correct Assembly Errors Using Linked Reads From Large Molecules"
author: [Shaun D Jackman^1^ (0000-0002-9275-5966), Justin Chu^1^, Rene L Warren^1^, Benjamin P Vandervalk^1^, Lauren Coombe^1^, Sarah Yeo^1^, Zhuyi Xue^1^, Hamid Mohamadi^1^, Joerg Bohlmann^2^, Steven JM Jones^1^, Inanc Birol^1^ (0000-0003-0950-7839)]
bibliography: tigmint.bib
csl: tigmint.csl
rangeDelim: "&ndash;"
eqnPrefix: "Equation"
figPrefix: "Fig."
tblPrefix: ["Table", "Tables"]
fontsize: 12pt
geometry: margin=1in
permalink: /
---

| ^1^ British Columbia Cancer Agency, Genome Sciences Centre, Vancouver, BC V5Z 4S6, Canada
| ^2^ University of British Columbia, Michael Smith Laboratories, Vancouver, BC V6T 1Z4, Canada

sjackman@bcgsc.ca; jchu@bcgsc.ca; rwarren@bcgsc.ca; benv@bcgsc.ca; lcoombe@bcgsc.ca; sarah.yeo@alumni.ubc.ca; zxue@bcgsc.ca; hmohamadi@bcgsc.ca; bohlmann@mail.ubc.ca; sjones@bcgsc.ca; ibirol@bcgsc.ca

> **Abstract.** Genome sequencing yields the sequence of many short snippets of DNA (reads) from a genome. Genome assembly attempts to reconstruct the original genome from which these reads were derived. This task is difficult due to gaps and errors in the sequencing data, repetitive sequence in the underlying genome, and heterozygosity, and assembly errors are common. These misassemblies may be identified by comparing the sequencing data to the assembly, and by looking for discrepancies between the two. Once identified, these misassemblies may be corrected, improving the quality of the assembly. Although tools exist to identify and correct misassemblies using Illumina pair-end and mate-pair sequencing, no such tool yet exists that makes use of the long distance information of the large molecules provided by linked reads, such as those offered by the 10x Genomics Chromium platform. We have developed the tool Tigmint for this purpose. To demonstrate the effectiveness of Tigmint, we corrected an assembly of the human genome using short reads assembled with ABySS 2.0. Tigmint reduced the number of breakpoints in the ABySS assembly by 319. While scaffolding with ARCS alone nearly doubled the NGA50 of the assembly at the scaffolding stage, the combination of Tigmint and ARCS tripled the scaffold NGA50 of the assembly to nearly 9 Mbp. This notable improvement in contiguity highlights the utility of assembly correction in refining assemblies. The source code of Tigmint is available for download from <https://github.com/bcgsc/tigmint>, and is distributed under the GNU GPL v3.0 license.

> **Keywords.** 10x Genomics Chromium reads &middot; De novo assembly &middot; Assembly correction &middot; Genome scaffolding &middot; Linked reads &middot; ABySS &middot; ARCS

# Introduction

Assemblies of short read sequencing data are easily confounded by repetitive sequences larger than the fragment size of the sequencing library. When the size of a repeat exceeds the library fragment size, the contig comes to an end in the best case, or results in misassembled sequence in the worst case. Misassemblies not only complicate downstream analyses, but also limit the contiguity of the assembly, when incorrectly assembled sequences prevent joining their adjacent and correctly assembled sequences during assembly scaffolding.

Long-read sequencing technologies have greatly improved assembly contiguity, by their ability to span these repeats, but at a cost roughly ten times that of short-read sequencing technology. For population studies and when sequencing large genomes, such as conifer genomes and other economically important crop species, this cost may be prohibitive. The 10x Genomics (Pleasanton, CA) Chromium technology generates linked reads from large DNA molecules at a cost comparable to standard short-read sequencing technologies. Whereas paired-end sequencing gives two reads from a small DNA fragment, linked reads yield roughly a hundred read pairs from molecules with a typical size of 10 to 100 kilobases. Linked reads indicate which reads were derived from the same DNA molecule (or molecules, when they share the same barcode), and so should be in close proximity in the underlying genome. The technology has been used previously to phase diploid genomes using a reference [@Zheng_2016], *de novo* assemble complex genomes in the gigabase scale [@Weisenfeld_2017], and further scaffold draft assemblies [@Mostovoy_2016].

A number of software tools employ linked reads for various applications. The LongRanger tool maps reads to repetitive sequence, phase small variants, and identify structural variants (<https://www.10xgenomics.com/software/>), and Supernova [@Weisenfeld_2017] assembles diploid genome sequences, both tools developed by the vendor. Among tools from academic labs, GROC-SVs [@Spies_2017], NAIBR [@Elyanow_2017], and Topsorter (<https://github.com/hanfang/Topsorter>) identify structural variants, and ARCS [@Yeo_2017], Architect [@Kuleshov_2016], and fragScaff [@Adey_2014] scaffold genome assemblies using linked reads.

In *de novo* sequencing projects, it is challenging yet important to ensure the correctness of the resulting assemblies. Tools to correct misassemblies typically inspect the reads aligned back to the assembly to identify discrepancies. Pilon [@Walker_2014] inspects the alignments to identify variants and correct small-scale misassemblies. NxRepair [@Murphy_2015] uses Illumina mate-pair sequencing to correct large-scale structural misassemblies. Linked reads offer an opportunity to use the long-range information provided by large molecules to identify misassemblies, yet no software tool currently exists to correct misassemblies using linked reads. Here we introduce a software tool, Tigmint, to identify misassemblies using this new and useful data type.

Tigmint first aligns reads to the assembly, and infers the extents of the large DNA molecules from these alignments. It then searches for atypical drops in physical molecule coverage. Since the physical coverage of the large molecules is more consistent and less prone to coverage dropouts than that of the short read sequencing data, it can be used to reveal the positions of possible misassemblies. Further, coincident clipped alignments of the first read of a molecule are used to refine the coordinates of detected misassemblies with base-pair accuracy. Finally, the linked reads may be used to scaffold the corrected assembly with ARCS [@Yeo_2017] and LINKS [@Warren_2015].

# Methods

## Algorithm

The data analysis pipeline of Tigmint is illustrated in @fig:pipeline and described below. The user provides a draft assembly in FASTA format and the reads in FASTQ format. Tigmint first aligns the reads to the reference using BWA-MEM [@Li_2013]. The alignments are filtered by mapping quality, alignment score, and number of mismatches to remove poorly aligned reads with the default thresholds $\textrm{MAPQ} > 0$, $\textrm{AS} \geq 100$, and $\textrm{NM} < 5$. Reads with the same barcode that map within 50,000 bp of the adjacent read are grouped into a molecule and assigned a unique numeric molecule identifier. A BED file is constructed, where each record indicates the start and end of one molecule, and records the number of reads that compose that molecule, their median mapping quality, alignment score, and number of mismatches. Unusually small molecules, shorter than 2000 bp by default, are filtered out.

![A flow chart of the data analysis pipeline of Tigmint. Inputs files are shown in parallelograms. Intermediate steps are shown in rectangles. Output files are shown in ellipses. File formats are indicated in parentheses. The optional steps of scaffolding and aligning to a reference genome to compute assembly metrics are shown in a dashed rectangle.](figures/pipeline.png){#fig:pipeline}

Regions with poor physical molecule coverage indicate potential misassemblies. The depth of molecule coverage at each position is computed from the molecule BED file using Bedtools [@Quinlan_2010]. The median and inter-quartile range (IQR) of the molecule depth of coverage is computed to determine the range of typical physical coverage for the experiment. Regions with coverage less than a threshold specified by the user, and suggested by default to be the median molecule coverage minus two times the IQR, are flagged as potential misassemblies.

The alignment to the assembly of the initial (left-most) read of a molecule that spans a misassembly will be clipped, and so the inferred genomic range of that molecule will start at precisely the position of the breakpoint. These clipped molecules are used to refine the breakpoint coordinates to base-pair resolution. Two or more molecules starting at the same position are required by default to identify a breakpoint, though this threshold parameter may be adjusted by the user.

The locations of these breakpoints are written to a tab-separated-values (TSV) file along with a summary of the evidence supporting the breakpoint, the number of molecules spanning that position (depth), and the number of molecules starting at that position (starts). The sequences of the original draft assembly are split at these breakpoints, producing a FASTA file. Three breakpoints in a region of half a megabase and the evidence used to identify them are illustrated in @fig:depth-starts.

![The depth of molecule coverage spanning each position, and the number of molecules whose first aligned nucleotide starts at each position. The horizontal dashed orange lines indicate the threshold parameters. The vertical dashed blue lines indicates the breakpoints identified by Tigmint.](figures/depth-starts.png){#fig:depth-starts}

Tigmint will optionally run ARCS [@Yeo_2017] and LINKS [@Warren_2015] at this point to scaffold these corrected sequences and improve the contiguity of the assembly. Tigmint will optionally align the scaffolds to a reference genome, if one is provided, to compute contiguity (NGA50) and correctness (number of breakpoints) of the assemblies before Tigmint, after Tigmint, and after ARCS and LINKS. The assembly metrics are calculated using ABySS-samtobreak, included with ABySS [@Jackman_2017]. Each breakpoint identified by ABySS-samtobreak indicates a difference between the assembly and the reference. These breakpoints are composed of both misassemblies and real structural variation between the reference and the sequenced genomes. We use NGA50, and its difference from NG50 as a surrogate for assembly correctness. The NGA50 metric summarizes both assembly contiguity and correctness by computing the NG50 of the lengths of alignment blocks to a reference genome, correcting the contiguity metric by accounting for potential misassemblies. It however also penalizes sequences at points of true variation between the sequenced and reference genomes. The true but unknown contiguity of the assembly, which accounts for misassemblies but not for structural variation, therefore lies somewhere between the lower bound of NGA50 and the upper bound of NG50.

## Human data set

We downloaded the ABySS 2.0 [@Jackman_2017] assembly `abyss-2.0/scaffolds.fa` from <http://bit.ly/ncbi-giab-abyss2> for the Genome in a Bottle (GIAB) HG004, assembled from Illumina paired-end and mate-pair reads [@Zook_2016]. We downloaded the 10x Genomics Chromium reads for this same individual from <http://bit.ly/giab-hg004-chromium> and used the LongRanger Basic pipeline to extract the barcodes from the reads. We ran Tigmint to correct the ABySS 2.0 assembly of HG004 using these Chromium reads with the parameters `depth_threshold=100 starts_threshold=2`. The choice of threshold parameters is discussed in the results. Both the uncorrected and corrected assembly are scaffolded using ARCS and LINKS. These assemblies are aligned to the GRCh38 reference genome using `bwa mem -xintractg`. The scaffold NGA50 and number of breakpoints are calculated using `abyss-samtobreak -G3088269832 -q10 -l500`.

This same data set was assembled by @Jackman_2017 using DISCOVARdenovo scaffolded using BESST [@Sahlin_2016]. We repeated this analysis of correction with Tigmint, scaffolding with ARCS, and comparison to the reference, using the DISCOVARdenovo + BESST assembly. We assembled these same linked reads of HG004 with Supernova 1.1.0 [@Weisenfeld_2017], and repeated this Tigmint + ARCS analysis using the Supernova assembly. The script to run this analysis is available online at <https://github.com/sjackman/tigmint-data>.

# Results

Correcting the ABySS assembly of the human data set HG004 with Tigmint reduces the number of breakpoints identified by ABySS-samtobreak by 319, a reduction of 7%. While the scaffold NG50 decreases slightly from 3.49 Mbp to 3.30 Mbp, the scaffold NGA50 remains unchanged; thus in this case, correcting the assembly with Tigmint improves the correctness of the assembly without substantially reducing its contiguity. However, scaffolding the uncorrected and corrected assemblies with ARCS yield markedly different results: nearly a two-fold increase (from 2.88 Mbp to 5.34 Mbp) versus a three-fold increase (to 8.69 Mbp) in NGA50, respectively. Further, correcting the assembly and then scaffolding yields a final assembly that is also more correct, as shown in @fig:metrics and @tbl:metrics.

Correcting the DISCOVARdenovo + BESST assembly reduces the number of breakpoints by 114 (a reduction of 2%). Using Tigmint to correct the assembly before scaffolding with ARCS yields an increase in NGA50 of 9% over using ARCS without Tigmint. As with the ABySS assembly, the DISCOVARdenovo + BESST + Tigmint + ARCS assembly is both more contiguous and has fewer breakpoints than the DISCOVARdenovo + BESST assembly. With two unrelated assembly pipelines, using Tigmint and ARCS together improves both the contiguity and correctness of the assembly.

Correcting the Supernova assembly of the HG004 linked reads with Tigmint reduces the number of breakpoints by 167 (a reduction of 3%), and after scaffolding the corrected assembly with ARCS, we see a slight reduction (2%) in NGA50 compared to the original Supernova assembly. Without Tigmint, the ABySS + ARCS assembly has a similar contiguity (scaffold NGA50) and fewer breakpoints than the Supernova assembly. The ABySS + Tigmint + ARCS assembly is both more contiguous and has fewer breakpoints compared to the Supernova assembly. Since the Supernova assembly is composed entirely of the linked reads, we do not expect significant gains from using these same data to correct the Supernova assembly. The Supernova assembly however has not made use of the mate-pair reads, and correcting the Supernova assembly with mate-pair reads may be an interesting area for future development of Tigmint.

![The assembly contiguity (scaffold NGA50) and correctness (number of breakpoints) metrics with and without correction using Tigmint prior to scaffolding with ARCS. The most contiguous and correct assemblies are found in the top-left corner. The DISCOVARdenovo + BESST assembly is labeled DISCOVAR.](figures/metrics.png){#fig:metrics}

Table: The assembly contiguity (scaffold NG50 and NGA50) and correctness (number of breakpoints) metrics with and without correction using Tigmint prior to scaffolding with ARCS. The reduction in the number of breakpoints from the row above it is shown in the final column. {#tbl:metrics}

| Assembly                          | NG50 (Mbp) | NGA50 (Mbp) | Breakpoints | Reduction |
| --------------------------------- | ---------: | ----------: | ----------: | --------: |
| ABySS                             |       3.49 |        2.88 |       4,770 |        NA |
| ABySS + Tigmint                   |       3.30 |        2.88 |       4,451 |       319 |
| ABySS + ARCS                      |       7.57 |        5.34 |       4,826 |        NA |
| ABySS + Tigmint + ARCS            |      11.54 |        8.69 |       4,485 |       341 |
| DISCOVAR + BESST                  |       6.92 |        3.57 |       5,680 |        NA |
| DISCOVAR + BESST + Tigmint        |       6.54 |        3.53 |       5,566 |       114 |
| DISCOVAR + BESST + ARCS           |      18.15 |        6.09 |       5,694 |        NA |
| DISCOVAR + BESST + Tigmint + ARCS |      19.39 |        6.62 |       5,577 |       117 |
| Supernova                         |      13.47 |        4.39 |       7,549 |        NA |
| Supernova + Tigmint               |       9.91 |        3.51 |       7,382 |       167 |
| Supernova + ARCS                  |      21.23 |        4.91 |       7,596 |        NA |
| Supernova + Tigmint + ARCS        |      14.20 |        4.31 |       7,437 |       159 |

The alignments of the ABySS assembly to the reference genome before and after Tigmint are visualized in @fig:jupiter using JupiterPlot (<https://github.com/JustinChu/JupiterPlot>), which makes use of Circos [@Krzywinski_2009]. The reference chromosomes are shown on the left in colour, and the assembly scaffolds are shown on the right in gray. The scaffolds on the right are arranged according the position of their best alignment to the reference. Chimeric scaffolds result in split alignments that manifest as lines criss-crossing the large coloured bands of concordant alignments. Small-scale structural variation is not visible due to the scale, but translocations (likely misassemblies) of sequences larger than 20 kbp are readily visible. A number of these split alignments are visible in the assembly before Tigmint, whereas after Tigmint no such split alignments are visible.

![The alignments to the reference genome of the ABySS assembly before and after Tigmint. Translocations are visible as lines criss-crossing the large coloured bands of concordant alignments. No translocations are visible after Tigmint.](figures/jupiter.png){#fig:jupiter}

The median molecule depth of this data is 163, and its inter-quartile range (IQR) is 31. We set the depth threshold parameter of Tigmint to 100, the median depth minus two times the IQR. The effect of varying the depth and starts threshold parameters of Tigmint on the assembly contiguity and correctness metrics is shown in @fig:parameters and @tbl:parameters. The assembly metrics are relatively insensitive to varying the depth threshold parameter. The starts threshold parameter specifies the number of molecules starting at the same position required to break the scaffold at that position. Multiple molecules starting at the same position is required to determine the position of the breakpoint. We tested thresholds of 2, 3, and 4 coincident molecules, and we observed the best performance when requiring two coincident molecules to break a scaffold.

![The effect of varying the depth and starts threshold parameters of Tigmint on the scaffold NGA50 and number of breakpoints. The original assembly is ABySS + ARCS without Tigmint. The most contiguous and correct assemblies are found in the top-left corner.](figures/parameters.png){#fig:parameters}

Table: The effect of varying the depth and starts threshold parameters of Tigmint on the scaffold NG50 and NGA50 and number of breakpoints. The first row is the original assembly, ABySS + ARCS without Tigmint. The reduction in the number of breakpoints from the original assembly is shown in the final column. {#tbl:parameters}

| Depth | Starts | NG50 (Mbp) | NGA50 (Mbp) | Breakpoints | Reduction |
| ----: | -----: | ---------: | ----------: | ----------: | --------: |
|    NA |     NA |       7.57 |        5.34 |       4,826 |        NA |
|   100 |      4 |       7.87 |        5.55 |       4,783 |        43 |
|   100 |      3 |       9.40 |        6.90 |       4,694 |       132 |
|    80 |      2 |      11.25 |        8.67 |       4,506 |       320 |
|   100 |      2 |      11.54 |        8.69 |       4,485 |       341 |
|   120 |      2 |      11.25 |        8.85 |       4,474 |       352 |

# Discussion

When aligning an assembly of an individual's genome to a reference genome of its species, we expect to see breakpoints where the assembled genome differs from the reference genome. These breakpoints are caused by both misassemblies and true differences between the individual and the reference. The median number of mobile-element insertions for example, just one class of structural variant, is estimated to be 1,218 per individual [@Sudmant_2015]. Misassemblies can be corrected by inspecting the alignments of the reads to the assembly. Correcting these misassemblies reduces the number of breakpoints when compared to the reference. Breakpoints due to true structural variation will however remain. For this reason, even a perfectly corrected assembly is expected to have a number of differences when compared to the reference.

Tigmint uses linked reads to reduce the number of misassemblies in a genome sequence assembly. The contiguity of the assembly is not appreciably affected by such a correction, while yielding an assembly that is more correct. Most scaffolding tools order and orient the sequences that they are given, but do not attempt to correct misassemblies. These misassemblies hold back the contiguity that can be achieved by scaffolding. Two sequences that should be connected together cannot be when one of those two sequences is connected incorrectly to a third sequence. By first correcting these misassemblies, the scaffolding tool can do a better job of connecting sequences, and we observe precisely this harmonious effect. Scaffolding an assembly that has been corrected with Tigmint yields a final assembly that is both more correct and substantially more contiguous than an assembly that has not been corrected.

# References
