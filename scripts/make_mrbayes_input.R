library(ape)

dir.create("results/trees/mrbayes", showWarnings = FALSE, recursive = TRUE)

dna <- read.dna("results/msa/M.aln.fa", format = "fasta")
chars <- as.character(dna)
seqs <- lapply(seq_len(nrow(chars)), function(i) chars[i, ])
names(seqs) <- rownames(chars)

out_file <- "results/trees/mrbayes/M.nex"

write.nexus.data(seqs, file = out_file, format = "dna", interleaved = FALSE)

mb_block <- c(
  "",
  "BEGIN MRBAYES;",
  "  set autoclose=yes nowarn=yes;",
  "  lset nst=6 rates=gamma;",
  "  prset statefreqpr=dirichlet(1,1,1,1);",
  "  mcmc ngen=100000 samplefreq=100 printfreq=100 diagnfreq=1000 nchains=4 nruns=2 savebrlens=yes filename=results/trees/mrbayes/M;",
  "  sump burnin=250;",
  "  sumt burnin=250;",
  "END;"
)

write(mb_block, file = out_file, append = TRUE)
