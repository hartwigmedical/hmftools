package com.hartwig.hmftools.sage.coverage;

import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import htsjdk.samtools.SAMRecord;

public class CoverageConsumer implements Consumer<SAMRecord> {

    private final String contig;
    private final List<GeneCoverage> geneCoverage;

    public CoverageConsumer(final String contig, final List<GeneCoverage> geneCoverage) {
        this.contig = contig;
        this.geneCoverage = geneCoverage.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
    }

    @Override
    public void accept(final SAMRecord record) {
        if (!record.getContig().equals(contig)) {
            throw new IllegalStateException();
        }

        final GenomeRegion alignment = GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
        geneCoverage.forEach(x -> x.accept(alignment));
    }
}
