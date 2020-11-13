package com.hartwig.hmftools.sage.coverage;

import java.util.Collection;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class GeneCoverage implements Consumer<GenomeRegion> {
    private static final int MAX_BUCKET = 37;

    private final String chromosome;
    private final String gene;
    private final List<ExonCoverage> exonCoverage;
    private final long minPosition;
    private final long maxPosition;

    public GeneCoverage(final List<NamedBed> exons) {
        assert (!exons.isEmpty());
        this.chromosome = exons.get(0).chromosome();
        this.gene = exons.get(0).name();
        this.exonCoverage = exons.stream().map(ExonCoverage::new).collect(Collectors.toList());

        long tmpMin = exons.get(0).start();
        long tmpMax = exons.get(0).end();

        for (NamedBed exon : exons) {
            tmpMin = Math.min(tmpMin, exon.start());
            tmpMax = Math.max(tmpMax, exon.end());
        }

        minPosition = tmpMin;
        maxPosition = tmpMax;
    }

    @Override
    public void accept(final GenomeRegion alignment) {
        if (alignment.chromosome().equals(chromosome)) {
            if (alignment.start() <= maxPosition && alignment.end() >= minPosition) {
                exonCoverage.forEach(x -> x.accept(alignment));
            }
        }
    }

    public String chromosome() {
        return chromosome;
    }

    @NotNull
    public GeneDepth geneDepth() {
        return ImmutableGeneDepth.builder().gene(gene).depthCounts(baseCoverageSummary(exonCoverage)).build();
    }

    static int[] baseCoverageSummary(Collection<ExonCoverage> exons) {
        int[] geneDepth = new int[MAX_BUCKET + 1];
        for (ExonCoverage exon : exons) {
            for (int baseDepth : exon.coverage()) {
                geneDepth[bucket(baseDepth)]++;
            }
        }

        return geneDepth;
    }

    static int bucket(int depth) {
        if (depth <= 30) {
            return depth;
        }

        for (int i = 0; i < 7; i++) {
            int maxDepth = 40 + 10 * i;
            if (depth < maxDepth) {
                return 30 + i;
            }

        }

        return MAX_BUCKET;
    }
}
