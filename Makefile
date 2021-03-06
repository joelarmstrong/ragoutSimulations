# Combine an evolver setup split into separate chromosomes into one
# large genome and get it ready for the evolverSimControl suite.
#
# Some things have been cribbed from Dent Earl's
# evolverInfileGeneration scripts.

# Settings.
rootName:=root
rootPath:=./root
chrs:=chr11 chr12 chr13 chr14 chr15 chr16 #chr17 chr18 chr19 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 
dataPath:=.

chrPaths:=$(foreach chr, $(chrs), ${dataPath}/${chr})
seqPaths:=$(foreach chr, $(chrs), ${dataPath}/${chr}/seq.rev)

%/DONE: %.nh
	simCtrl_runSim.py --inputNewick $(shell cat $<) --outDir $(dir $@) --rootDir root --paramsDir params --maxThreads 30 --maxCpu 30 --noMEs --jobTree jobTree$(shell basename $<)
	touch $@

all: ${rootPath}/annots.gff ${rootPath}/seq.rev stats

${rootPath}/annots.gff: $(foreach chrPath, ${chrPaths}, ${chrPath}/annots.gff.new)
	@mkdir -p $(dir $@)
# evolver_evo breaks when parsing the GFF lines that don't have a
# gene_index, so we get rid of any entries that don't have an
# index. This unfortunately means the non-coding conserved elements
# are currently not included.
	grep -h gene_index $^ > $@.tmp
	mv $@.tmp $@

${rootPath}/seq.rev: ${rootPath}/root.fa
	@mkdir -p $(dir $@)
	evolver_cvt -fromfasta $< -genome ${rootName} -torev $@.tmp
	mv $@.tmp $@

${rootPath}/root.fa: $(foreach chrPath, ${chrPaths}, ${chrPath}/seq.fa)
	@mkdir -p $(dir $@)
	cat $^ > $@.tmp
	mv $@.tmp $@

$(foreach chr, $(chrs), ${dataPath}/${chr}/seq.fa): %/seq.fa: %/seq.rev
	@mkdir -p $(dir $@)
	evolver_cvt -fromrev $< -tofasta $@.tmp
	mv $@.tmp $@

# The following is taken pretty much verbatim from evolverInfileGeneration.

stats: ${rootPath}/stats/expanded_annots.gff ${rootPath}/stats/annotstats.txt

${rootPath}/stats/annotstats.txt: ${rootPath}/annots.gff ${rootPath}/seq.rev
	@mkdir -p $(dir $@)
	evolver_evo -annotstats ${rootPath}/annots.gff -seq ${rootPath}/seq.rev -log $@.tmp
	mv $@.tmp $@

${rootPath}/stats/expanded_annots.gff: ${rootPath}/stats/introns.gff ${rootPath}/stats/exons.gff ${rootPath}/annots.gff
	@mkdir -p $(dir $@)
	cat $^ > $@.tmp
	mv $@.tmp $@

${rootPath}/stats/introns.gff: ${rootPath}/stats/exons.gff ${rootPath}/stats/cds_annots.gff ${rootPath}/annots.gff
	@mkdir -p $(dir $@)
	evolver_gff_exons2introns.py $< > $@.tmp
	mv $@.tmp $@

${rootPath}/stats/exons.gff: ${rootPath}/stats/cds_annots.gff ${rootPath}/annots.gff
	@mkdir -p $(dir $@)
	evolver_gff_cdsutr2exons.py $< > $@.tmp
	mv $@.tmp $@

${rootPath}/stats/cds_annots.gff: ${rootPath}/annots.gff ${rootPath}/stats/merged_root.stats.txt ${rootPath}/stats/merged_branch.stats.txt
	@mkdir -p $(dir $@)
	(egrep 'CDS|UTR' $< || true) | cat > $@.tmp
	mv $@.tmp $@

${rootPath}/stats/merged_root.stats.txt:
	@mkdir -p $(dir $@)
	touch $@

${rootPath}/stats/merged_branch.stats.txt:
	@mkdir -p $(dir $@)
	touch $@
