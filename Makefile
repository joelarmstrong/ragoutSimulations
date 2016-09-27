# Combine an evolver setup split into separate chromosomes into one
# large genome and get it ready for the evolverSimControl suite.
#
# Some things have been cribbed from Dent Earl's
# evolverInfileGeneration scripts.

# Settings.
rootName:=root
rootPath:=./root
chrs:=chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
dataPath:=.

chrPaths:=$(foreach chr, $(chrs), ${dataPath}/${chr})
seqPaths:=$(foreach chr, $(chrs), ${dataPath}/${chr}/seq.rev)

all: ${rootPath}/annots.gff ${rootPath}/seq.rev

${rootPath}/annots.gff: $(foreach chrPath, ${chrPaths}, ${chrPath}/annots.gff.new)
	@mkdir -p $(dir $@)
	cat $^ > $@.tmp
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
