package com.hartwig.hmftools.common.purple;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosomes;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.cnv.ImmutableCNVCopyNumber;

import java.util.List;
import java.util.function.Consumer;

public enum PadCopyNumber {
    ;

    public static List<CopyNumber> pad(List<CopyNumber> copyNumbers) {
        PadCopyNumberForChromosome padder = null;
        List<CopyNumber> result = Lists.newArrayList();
        for (CopyNumber copyNumber : copyNumbers) {
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

    static class PadCopyNumberForChromosome {

        private final String chromosome;
        private final Consumer<CopyNumber> handler;
        private long currentPosition = 1;

        PadCopyNumberForChromosome(String chromosome, Consumer<CopyNumber> handler) {
            this.chromosome = chromosome;
            this.handler = handler;
        }

        String getChromosome() {
            return chromosome;
        }

        void addCopyNumber(CopyNumber number) {
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

        private CopyNumber createDefaultCopyNumber(long aStart, long aEnd) {
            return ImmutableCNVCopyNumber.of(2, chromosome, aStart, aEnd, null);
        }
    }
}
