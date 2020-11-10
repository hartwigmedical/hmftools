package com.hartwig.hmftools.sage.coverage;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

class ExonCoverage implements GenomeRegion, Consumer<GenomeRegion> {

    private final NamedBed exon;
    private final int[] coverage;

    public ExonCoverage(final NamedBed exon) {
        this.exon = exon;
        this.coverage = new int[(int) exon.bases()];
    }

    @Override
    public void accept(final GenomeRegion alignment) {
        if (!alignment.chromosome().equals(exon.chromosome())) {
            return;
        }

        synchronized (coverage) {
            int start = (int) Math.max(start(), alignment.start());
            int end = (int) Math.min(end(), alignment.end());

            for (int i = index(start); i <= index(end); i++) {
                coverage[i] += 1;
            }
        }
    }

    @NotNull
    public int[] coverage() {
        return coverage;
    }

    @NotNull
    public String gene() {
        return exon.name();
    }

    @NotNull
    @Override
    public String chromosome() {
        return exon.chromosome();
    }

    @Override
    public long start() {
        return exon.start();
    }

    @Override
    public long end() {
        return exon.end();
    }

    private int index(int position) {
        return (int) (position - start());
    }
}
