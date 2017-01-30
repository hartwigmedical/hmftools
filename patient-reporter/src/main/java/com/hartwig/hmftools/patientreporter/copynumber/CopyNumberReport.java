package com.hartwig.hmftools.patientreporter.copynumber;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CopyNumberReport {

    @VisibleForTesting
    static final String COPY_NUMBER_GAIN = "gain";
    @VisibleForTesting
    static final String COPY_NUMBER_LOSS = "loss";
    @VisibleForTesting
    static final String COPY_NUMBER_NEUTRAL = "none";

    @NotNull
    private final String gene;
    @NotNull
    private final String transcript;
    private final int copyNumber;

    private CopyNumberReport(@NotNull final String gene, @NotNull final String transcript, final int copyNumber) {
        this.gene = gene;
        this.transcript = transcript;
        this.copyNumber = copyNumber;
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

    public static class Builder {
        @NotNull
        private String gene = Strings.EMPTY;
        @NotNull
        private String transcript = Strings.EMPTY;
        private int copyNumber = 0;

        public Builder() {
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
            return new CopyNumberReport(gene, transcript, copyNumber);
        }
    }
}
