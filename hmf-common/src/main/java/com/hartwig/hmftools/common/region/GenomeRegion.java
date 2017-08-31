package com.hartwig.hmftools.common.region;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public interface GenomeRegion extends Comparable<GenomeRegion> {

    @NotNull
    String chromosome();

    long start();

    long end();

    default long bases() {
        return 1 + end() - start();
    }

    @Override
    default int compareTo(@NotNull final GenomeRegion other) {

        if (chromosome().equals(other.chromosome())) {
            if (start() < other.start()) {
                return -1;
            } else if (start() == other.start()) {
                return 0;
            }
            return 1;
        }

        return HumanChromosome.fromString(chromosome()).compareTo(HumanChromosome.fromString(other.chromosome()));
    }

    default boolean overlaps(@NotNull final GenomeRegion other) {
        return other.chromosome().equals(chromosome()) && other.end() > start() && other.start() < end();
        // TODO: Double check this could be replaced with "return overlappingBases(other) > 0;"
    }

    default long overlappingBases(@NotNull final GenomeRegion other) {
        if (!chromosome().equals(other.chromosome())) {
            return 0;
        }

        long minEnd = Math.min(end(), other.end());
        long maxStart = Math.max(start(), other.start());
        return Math.max(0, 1 + minEnd - maxStart);
    }
}
