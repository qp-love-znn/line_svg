# line_svg
get the graph accroding to the position of circ and gene 
run it as follows : perl graph_main.pl line_svg -conf circ_gene.conf -type circ_gene -output out
the circ_gene.conf like this :
    width = 1600
    height = 1600
    cir_posit = inputdir/circ_count_bigTo5.bed
    gene_posit = inputdir/annot.txt
    chr_length = 15222
   
the circ_count_bigTo5.bed format like this :
  circ_ID\tstart\tend
the circ_gene format like this :
  gene_ID\tstart\tend
