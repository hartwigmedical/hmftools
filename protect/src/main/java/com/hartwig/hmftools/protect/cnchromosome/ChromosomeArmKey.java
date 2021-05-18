package com.hartwig.hmftools.protect.cnchromosome;

import java.util.Objects;

import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.jetbrains.annotations.NotNull;

public class ChromosomeArmKey {

    @NotNull
    private final String chromosome;
    @NotNull
    private final ChromosomeArm chromosomeArm;

    public ChromosomeArmKey(@NotNull final String chromosome, @NotNull final ChromosomeArm chromosomeArm) {
        this.chromosome = chromosome;
        this.chromosomeArm = chromosomeArm;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final ChromosomeArmKey that = (ChromosomeArmKey) o;
        return chromosome.equals(that.chromosome) && chromosomeArm == that.chromosomeArm;
    }

    @Override
    public int hashCode() {
        return Objects.hash(chromosome, chromosomeArm);
    }

    @Override
    public String toString() {
        return "ChromosomeArmKey{" + "chromosome='" + chromosome + '\'' + ", chromosomeArm=" + chromosomeArm + '}';
    }
}
