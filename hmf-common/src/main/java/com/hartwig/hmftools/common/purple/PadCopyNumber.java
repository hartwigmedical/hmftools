package com.hartwig.hmftools.common.purple;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.copynumber.freec.FreecCopyNumber;
import com.hartwig.hmftools.common.copynumber.freec.ImmutableFreecCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public enum PadCopyNumber {
    ;

    @NotNull
    public static List<GenomeRegion> pad(@NotNull final List<GenomeRegion> copyNumbers) {
        PadCopyNumberForChromosome padder = null;
        final List<GenomeRegion> result = Lists.newArrayList();
        for (final GenomeRegion copyNumber : copyNumbers) {
            if (padder == null) {
                padder = new PadCopyNumberForChromosome(copyNumber.chromosome(), result::add);
            } else if (!padder.chromosome.equals(copyNumber.chromosome())) {
                padder.complete();
                padder = new PadCopyNumberForChromosome(copyNumber.chromosome(), result::add);
            }
            padder.addCopyNumber(copyNumber);
        }

        return result;
    }

    private static class PadCopyNumberForChromosome {

        @NotNull
        private final String chromosome;
        @NotNull
        private final Consumer<GenomeRegion> handler;
        private long currentPosition = 1;

        PadCopyNumberForChromosome(@NotNull final String chromosome, @NotNull final Consumer<GenomeRegion> handler) {
            this.chromosome = chromosome;
            this.handler = handler;
        }

        @NotNull
        String getChromosome() {
            return chromosome;
        }

        void addCopyNumber(@NotNull final GenomeRegion number) {
            if (number.start() > currentPosition) {
                handler.accept(createDefaultCopyNumber(currentPosition, number.start() - 1));
            }

            handler.accept(number);
            currentPosition = number.end() + 1;
        }

        void complete() {
            long chromosomeLength = Chromosomes.length(chromosome);
            if (currentPosition < chromosomeLength) {
                handler.accept(createDefaultCopyNumber(currentPosition, chromosomeLength));
            }
        }

        @NotNull
        private FreecCopyNumber createDefaultCopyNumber(long start, long end) {
            return ImmutableFreecCopyNumber.builder().chromosome(chromosome).start(start).end(end).value(2).build();
        }
    }
}
