#test CC genoprobs
tr.id = "ENSMUSG00000028619"
marker.chr = 4
marker.pos = 107159244
nearest.geno <- pull_genoprobpos(cc.geno, map = cc.map, chr = marker.chr, pos = marker.pos)[u_strain_idx,]

#CC051xCC059 and CC012xCC002 should both be genotype EE

