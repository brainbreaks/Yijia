library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(GenomicRanges)

macs_cols = cols(
  macs_chrom=col_character(), macs_start=col_double(), macs_end=col_double(), macs_length=col_character(), macs_summit_abs=col_double(),
  macs_pileup=col_double(), macs_pvalue=col_double(), macs_fc=col_double(), macs_qvalue=col_double(), macs_name=col_character(), macs_comment=col_character()
)

tlx_cols = cols(
  Qname=readr::col_character(), JuncID=readr::col_character(), Rname=readr::col_character(), Junction=readr::col_double(),
  Strand=readr::col_character(), Rstart=readr::col_double(), Rend=readr::col_double(),
  B_Rname=readr::col_character(), B_Rstart=readr::col_double(), B_Rend=readr::col_double(), B_Strand=readr::col_double(),
  B_Qstart=readr::col_double(), B_Qend=readr::col_double(), Qstart=readr::col_double(), Qend=readr::col_double(), Qlen=readr::col_double(),
  B_Cigar=readr::col_character(), Cigar=readr::col_character(), Seq=readr::col_character(), J_Seq=readr::col_character(), Barcode=readr::col_logical(),
  unaligned=readr::col_double(), baitonly=readr::col_double(), uncut=readr::col_double(), misprimed=readr::col_double(), freqcut=readr::col_double(),
  largegap=readr::col_double(), mapqual=readr::col_double(), breaksite=readr::col_double(), sequential=readr::col_double(), repeatseq=readr::col_double(), duplicate=readr::col_double()
)

repeatmasker_cols = cols(
  repeatmasker_score=col_double(), repeatmasker_chrom=col_character(), repeatmasker_start=col_double(),
  repeatmasker_end=col_double(), repeatmasker_strand=col_character(), repeatmasker_name=col_character(),
  repeatmasker_class=col_character(), repeatmasker_family=col_character()
)


macs2 = function(sample, control, qvalue=0.01, extsize=200, output_dir="data/macs2", llocal=10000000) {
  bed_sample = sample
  bed_control = "data/JJ02_B400_012.bed"
  name = gsub(".bed$", "", basename(sample))
  cmd = stringr::str_glue("macs2 callpeak -t {bed_sample} -c {bed_control} --seed 123 -f BED -g hs --keep-dup all -n {name} --outdir {output_dir} --nomodel --extsize {extsize} -q {qvalue} --llocal {llocal} --bdg --trackline", bed_sample=sample, bed_control=control, name=name, output_dir=output_dir, extsize=extsize, qvalue=qvalue, llocal=sprintf("%0.0f", llocal))
  print(cmd)
  system(cmd)
}

main = function() {
  #TEst
  bait_region = 1e6
  baits_df = readr::read_tsv("data/baits.tsv")
  samples_df = readr::read_tsv("data/samples.tsv")
  repeatmasker_df = readr::read_tsv("data/ucsc_hg19_repeatmasker.tsv", col_names=names(repeatmasker_cols$cols), col_types=repeatmasker_cols, skip=1) %>%
    dplyr::mutate(repeatmasker_id=1:n())
  repeatmasker_ranges = GenomicRanges::makeGRangesFromDataFrame(repeatmasker_df %>% dplyr::mutate(seqnames=repeatmasker_chrom, start=repeatmasker_start, end=repeatmasker_end), keep.extra.columns=T)

  # system("singularity pull docker://sandrejev/htgts:latest")
  # system("singularity exec -B `pwd` htgts_latest.sif download hg19")

  tlx2repeatmasker_all = data.frame()
  for(f in list.files("data/htgts", pattern="*._result.tlx", full.names=T))
  {
    sample_file = gsub("_result.tlx", "", basename(f))

    tlx_df = readr::read_tsv(f, comment="#", skip=16, col_names=names(tlx_cols$cols), col_types=tlx_cols) %>%
      dplyr::mutate(tlx_id=1:n())

    tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Rstart, end=Rend), keep.extra.columns=T, ignore.strand=T)
    tlx2repeatmasker = as.data.frame(IRanges::findOverlaps(tlx_ranges, repeatmasker_ranges)) %>%
      dplyr::inner_join(repeatmasker_df, by=c("subjectHits"="repeatmasker_id")) %>%
      dplyr::group_by(queryHits) %>%
      dplyr::summarise(repeatmasker_name=paste(unique(repeatmasker_name), collapse=", "), repeatmasker_class=paste(unique(repeatmasker_class), collapse=", "), repeatmasker_family=paste(unique(repeatmasker_family), collapse=", "), repeatmasker_length=max(abs(repeatmasker_end-repeatmasker_start))) %>%
      dplyr::right_join(tlx_df, by=c("queryHits"="tlx_id")) %>%
      dplyr::mutate(sample_file=sample_file) %>%
      dplyr::inner_join(samples_df, by="sample_file") %>%
      dplyr::inner_join(baits_df, by="bait_id")
    tlx2repeatmasker_all = dplyr::bind_rows(tlx2repeatmasker_all, tlx2repeatmasker)

    f_norepeats = paste0(dirname(f), "/tmp/", gsub("_result.tlx", "_result_clean.tlx", basename(f)))
    f_bed = paste0(dirname(f), "/tmp/", gsub("_result.tlx", "_result_clean.bed", basename(f)))
    dir.create(file.path(dirname(f), "tmp"), showWarnings=F)

    tlx2repeatmasker_filter = tlx2repeatmasker %>%
      dplyr::filter(is.na(repeatmasker_class)) %>%
      dplyr::filter(bait_chrom==Rname & abs(Junction-bait_start)>=bait_region/2)
    readr::write_tsv(tlx2repeatmasker_filter[,names(tlx_cols$cols)], file=f_norepeats, na="")
    system(stringr::str_glue("singularity exec -B `pwd` htgts_latest.sif tlx2BED-MACS.pl {tlx} {bed} 0", tlx=f_norepeats, bed=f_bed))

    # system("singularity exec -B `pwd` htgts_latest.sif TranslocPreprocess.pl tutorial_metadata.txt preprocess --read1 pooled_R1.fastq.gz --read2 pooled_R2.fastq.gz")
    # system("singularity exec -B `pwd` htgts_latest.sif TranslocWrapper.pl tutorial_metadata.txt preprocess/ results/ --threads 2")
    # system(stringr::str_glue("macs2 pileup -i {bed} -f BED -o {output_name} --outdir {output_dir} --extsize 2000", bed=f_bed, output_name=sample, output_dir="data/macs2"))
    # system(stringr::str_glue("macs2 callpeak -t {bed} -f BED -g hs --keep-dup all -n {output_name} --outdir {output_dir} --nomodel --extsize 2000 -q 0.001 --llocal 10000000", bed=f_bed, output_name=sample, output_dir="data/macs2"))
    # macs2 callpeak -t {input} -f BED -g mm --keep-dup all -n {output} --outdir {outdir} --nomodel --extsize 2000 -q 0.01 --llocal 10000000
  }

  tlx2repeatmasker_ggplot = dplyr::bind_rows(
    tlx2repeatmasker_all %>% dplyr::mutate(filter="All"),
    tlx2repeatmasker_all %>% dplyr::mutate(filter="No repeats") %>% dplyr::filter(is.na(repeatmasker_class)),
    tlx2repeatmasker_all %>% dplyr::mutate(filter="All (excluding bait)") %>% dplyr::filter(bait_chrom==Rname & abs(Junction-bait_start)>=bait_region/2),
    tlx2repeatmasker_all %>% dplyr::mutate(filter="No repeats (excluding bait)") %>% dplyr::filter(is.na(repeatmasker_class) & bait_chrom==Rname & abs(Junction-bait_start)>=bait_region/2)
  ) %>% dplyr::filter(Rname=="chr1")
  samples_colors = with(tlx2repeatmasker_ggplot %>% dplyr::distinct(sample_desc, .keep_all=T), setNames(as.character(sample_color), sample_desc))

  tlx2repeatmasker_ggplot %>%
    reshape2::dcast(filter ~ sample_desc) %>%
    readr::write_tsv("meeting/reads_counts.tsv")

  tlx2repeatmasker_ggplot %>%
    dplyr::filter(filter=="All (excluding bait)") %>%
    tidyr::separate_rows(repeatmasker_class) %>%
    reshape2::dcast(repeatmasker_class ~ sample_desc) %>%
    readr::write_tsv("meeting/reads_counts_repeats.tsv")

  # LTR - identical sequence at end of retrotransposon
  # Satellite - tandem repeats
  # LINE 7000bp transposons


  ggplot(tlx2repeatmasker_ggplot) +
    geom_density(aes(x=Junction, color=sample_desc), bw=1e5, alpha=0.2) +
    geom_vline(aes(xintercept=bait_start-5e5), data=tlx2repeatmasker_ggplot %>% dplyr::distinct(bait_start)) +
    geom_vline(aes(xintercept=bait_end+5e5), data=tlx2repeatmasker_ggplot %>% dplyr::distinct(bait_end)) +
    scale_color_manual(values=samples_colors) +
    facet_wrap(~filter, ncol=2) +
    theme_bw(base_size=16)

  #
  # macs2(sample="data/JJ02_B400_012.bed", control="data/JJ01_B400_012.bed")
  # macs2(sample="data/JJ03_B400_012.bed", control="data/JJ01_B400_012.bed")
  macs2(sample="data/htgts/tmp/JJ03_B400_012_result_clean.bed", control="data/htgts/tmp/JJ01_B400_012_result_clean.bed", qvalue=0.01, llocal=1e6)
  macs2(sample="data/htgts/tmp/JJ02_B400_012_result_clean.bed", control="data/htgts/tmp/JJ01_B400_012_result_clean.bed", qvalue=0.01, llocal=1e6)
}



