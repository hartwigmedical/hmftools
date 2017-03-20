package com.hartwig.hmftools.patientreporter.copynumber;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CopyNumberReport implements Comparable<CopyNumberReport> {

    @VisibleForTesting
    static final String COPY_NUMBER_GAIN = "copy-gain";
    @VisibleForTesting
    static final String COPY_NUMBER_LOSS = "copy-loss";
    @VisibleForTesting
    static final String COPY_NUMBER_NEUTRAL = "none";

    @NotNull
    private final String chromosome;
    @NotNull
    private final String gene;
    @NotNull
    private final String transcript;
    private final int copyNumber;

    private CopyNumberReport(@NotNull final String chromosome, @NotNull final String gene,
            @NotNull final String transcript, final int copyNumber) {
        this.chromosome = chromosome;
        this.gene = gene;
        this.transcript = transcript;
        this.copyNumber = copyNumber;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @NotNull
    public String gene() {
        return gene;
    }

    @NotNull
    public String transcript() {
        return transcript;
    }

    public int copyNumber() {
        return copyNumber;
    }

    @NotNull
    public String resolveType() {
        if (copyNumber > 2) {
            return COPY_NUMBER_GAIN;
        } else if (copyNumber < 2) {
            return COPY_NUMBER_LOSS;
        } else {
            return COPY_NUMBER_NEUTRAL;
        }
    }

    @Override
    public int compareTo(@NotNull final CopyNumberReport other) {
        final Integer intChrom1 = toInteger(chromosome);
        final Integer intChrom2 = toInteger(other.chromosome);
        if (intChrom1 == null && intChrom2 == null) {
            return chromosome.compareTo(other.chromosome);
        } else if (intChrom1 == null) {
            return 1;
        } else if (intChrom2 == null) {
            return -1;
        } else {
            return intChrom1.compareTo(intChrom2);
        }
    }

    @Nullable
    private static Integer toInteger(@NotNull String string) {
        try {
            return Integer.parseInt(string);
        } catch (NumberFormatException exception) {
            return null;
        }
    }

    public static class Builder {
        @NotNull
        private String chromosome = Strings.EMPTY;
        @NotNull
        private String gene = Strings.EMPTY;
        @NotNull
        private String transcript = Strings.EMPTY;
        private int copyNumber = 0;

        public Builder() {
        }

        @NotNull
        public Builder chromosome(@NotNull final String chromosome) {
            this.chromosome = chromosome;
            return this;
        }

        @NotNull
        public Builder gene(@NotNull final String gene) {
            this.gene = gene;
            return this;
        }

        @NotNull
        public Builder transcript(@NotNull final String transcript) {
            this.transcript = transcript;
            return this;
        }

        public Builder copyNumber(final int copyNumber) {
            this.copyNumber = copyNumber;
            return this;
        }

        @NotNull
        public CopyNumberReport build() {
            return new CopyNumberReport(chromosome, gene, transcript, copyNumber);
        }
    }
}
