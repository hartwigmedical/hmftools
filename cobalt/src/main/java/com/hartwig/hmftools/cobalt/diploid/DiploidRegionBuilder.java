package com.hartwig.hmftools.cobalt.diploid;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

class DiploidRegionBuilder implements Consumer<DiploidCount> {

    private static final int WINDOW_SIZE = 1000;

    private final List<GenomeRegion> result = Lists.newArrayList();

    private final double cutoff;
    private final int maleSamples;
    private final int femaleSamples;
    private long totalDiploidBases = 0;

    private String contig = "";
    private long start = 0;
    private long end = 0;
    private boolean isDiploid;

    DiploidRegionBuilder(final double cutoff, final int femaleSamples, final int maleSamples) {
        this.cutoff = cutoff;
        this.maleSamples = maleSamples;
        this.femaleSamples = femaleSamples;
    }

    @Override
    public void accept(@NotNull DiploidCount count) {
        int samples = count.chromosome().equals("Y") || count.chromosome().equals("chrY") ? maleSamples : femaleSamples;
        boolean isCountDiploid = Doubles.greaterOrEqual(count.proportionIsDiploid(samples), cutoff);
        if (!count.chromosome().equals(contig) || isCountDiploid != isDiploid) {
            finaliseCurrent();
            contig = count.chromosome();
            start = count.position();
            isDiploid = isCountDiploid;
        }

        end = count.position() + WINDOW_SIZE - 1;
    }

    private void finaliseCurrent() {
        if (isDiploid) {
            result.add(GenomeRegions.create(contig, start, end));
            totalDiploidBases += end - start + 1;
        }
        isDiploid = false;
    }

    public long getTotalDiploidBases() {
        return totalDiploidBases;
    }

    @NotNull
    public List<GenomeRegion> build() {
        finaliseCurrent();
        return result;
    }
}
