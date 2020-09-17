package com.hartwig.hmftools.common.genome.chromosome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class CobaltChromosomes {

    static final int MIN_RATIO_COUNT = 1000;
    private static final double TWO_X_CUTOFF = 0.65;
    private static final double Y_CUTOFF = 0.05;
    private static final double MOSIAC_CUTOFF = 0.2;

    private final List<Chromosome> chromosomeList;
    private final Map<String, CobaltChromosome> chromosomeMap;
    private final boolean isFemale;

    public CobaltChromosomes(final List<MedianRatio> ratios) {
        this.chromosomeList = Lists.newArrayList();
        this.chromosomeMap = Maps.newHashMap();

        final double yMedian = contigRatio("Y", ratios);
        final double xMedian = contigRatio("X", ratios);
        isFemale = Doubles.greaterOrEqual(xMedian, TWO_X_CUTOFF) && Doubles.lessOrEqual(yMedian, Y_CUTOFF);

        for (MedianRatio ratio : ratios) {
            final String contig = ratio.chromosome();
            boolean isX = contig.equals("X") || contig.equals("chrX");
            boolean isY = contig.equals("Y") || contig.equals("chrY");

            int typicalCopies = typical(ratio.chromosome());

            double ratioImpliedCopies = ratio.medianRatio() * 2;
            double distanceFromCommonRatio = Math.abs(ratioImpliedCopies - Math.round(ratioImpliedCopies)) / 2.0;
            boolean isMosiac = Doubles.greaterThan(distanceFromCommonRatio, MOSIAC_CUTOFF);
            if (typicalCopies > 0) {
                CobaltChromosome chromosome = ImmutableCobaltChromosome.builder()
                        .contig(ratio.chromosome())
                        .impliedCopies((int) Math.round(ratioImpliedCopies))
                        .isAllosome(isX || isY)
                        .isAutosome(!isX && !isY)
                        .mosiac(isMosiac)
                        .build();

                chromosomeMap.put(contig, chromosome);
            }
        }
    }

    public int typical(String contig) {
        boolean isX = contig.equals("X") || contig.equals("chrX");
        boolean isY = contig.equals("Y") || contig.equals("chrY");

        if (isMale() && (isX || isY)) {
            return 1;
        }

        if (isFemale() && isY) {
            return 0;
        }

        return 2;
    }

    public boolean isFemale() {
        return isFemale;
    }

    public boolean isMale() {
        return !isFemale;
    }

    public boolean contains(@NotNull final String contig) {
        return chromosomeMap.containsKey(contig);
    }

    @NotNull
    public CobaltChromosome get(@NotNull final String contig) {
        return chromosomeMap.get(contig);
    }

    @NotNull
    public List<Chromosome> chromosomes() {
        return chromosomeList;
    }

    private double contigRatio(@NotNull final String contig, @NotNull final List<MedianRatio> ratios) {
        return ratios.stream()
                .filter(x -> x.count() >= MIN_RATIO_COUNT)
                .filter(x -> x.chromosome().equals(contig) || x.chromosome().equals("chr" + contig))
                .mapToDouble(MedianRatio::medianRatio)
                .findFirst()
                .orElse(0);
    }

}
