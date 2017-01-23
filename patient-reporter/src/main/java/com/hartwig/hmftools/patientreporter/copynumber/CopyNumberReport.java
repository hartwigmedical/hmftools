package com.hartwig.hmftools.patientreporter.copynumber;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CopyNumberReport {
    @NotNull
    private final String gene;
    @NotNull
    private final String transcript;
    @NotNull
    private final String finding;

    private CopyNumberReport(@NotNull final String gene, @NotNull final String transcript,
            @NotNull final String finding) {
        this.gene = gene;
        this.transcript = transcript;
        this.finding = finding;
    }

    @NotNull
    public String getGene() {
        return gene;
    }

    @NotNull
    public String getTranscript() {
        return transcript;
    }

    @NotNull
    public String getFinding() {
        return finding;
    }

    static class Builder {
        @NotNull
        private String gene = Strings.EMPTY;
        @NotNull
        private String transcript = Strings.EMPTY;
        @NotNull
        private String finding = Strings.EMPTY;

        Builder() {
        }

        @NotNull
        Builder gene(@NotNull final String gene) {
            this.gene = gene;
            return this;
        }

        @NotNull
        Builder transcript(@NotNull final String transcript) {
            this.transcript = transcript;
            return this;
        }

        @NotNull
        Builder finding(@NotNull final String finding) {
            this.finding = finding;
            return this;
        }

        @NotNull
        CopyNumberReport build() {
            return new CopyNumberReport(gene, transcript, finding);
        }
    }
}
