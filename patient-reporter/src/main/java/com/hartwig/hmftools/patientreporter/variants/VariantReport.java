package com.hartwig.hmftools.patientreporter.variants;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class VariantReport {

    @NotNull
    private final String gene;
    @NotNull
    private final String position;
    @NotNull
    private final String ref;
    @NotNull
    private final String alt;
    @NotNull
    private final String transcript;
    @NotNull
    private final String hgvsCoding;
    @NotNull
    private final String hgvsProtein;
    @NotNull
    private final String consequence;
    @NotNull
    private final String cosmicID;
    private final int totalReadCount;
    private final int alleleReadCount;

    private VariantReport(@NotNull final String gene, @NotNull final String position, @NotNull final String ref,
            @NotNull final String alt, @NotNull final String transcript, @NotNull final String hgvsCoding,
            @NotNull final String hgvsProtein, @NotNull final String consequence, @NotNull final String cosmicID,
            final int totalReadCount, final int alleleReadCount) {
        this.gene = gene;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.transcript = transcript;
        this.hgvsCoding = hgvsCoding;
        this.hgvsProtein = hgvsProtein;
        this.consequence = consequence;
        this.cosmicID = cosmicID;
        this.totalReadCount = totalReadCount;
        this.alleleReadCount = alleleReadCount;
    }

    @NotNull
    public String gene() {
        return gene;
    }

    @NotNull
    public String position() {
        return position;
    }

    @NotNull
    public String ref() {
        return ref;
    }

    @NotNull
    public String alt() {
        return alt;
    }

    @NotNull
    public String transcript() {
        return transcript;
    }

    @NotNull
    public String hgvsCoding() {
        return hgvsCoding;
    }

    @NotNull
    public String hgvsProtein() {
        return hgvsProtein;
    }

    @NotNull
    public String consequence() {
        return consequence;
    }

    @NotNull
    public String cosmicID() {
        return cosmicID;
    }

    public int totalReadCount() {
        return totalReadCount;
    }

    public int alleleReadCount() {
        return alleleReadCount;
    }

    public static class Builder {
        @NotNull
        private String gene = Strings.EMPTY;
        @NotNull
        private String position = Strings.EMPTY;
        @NotNull
        private String ref = Strings.EMPTY;
        @NotNull
        private String alt = Strings.EMPTY;
        @NotNull
        private String transcript = Strings.EMPTY;
        @NotNull
        private String hgvsCoding = Strings.EMPTY;
        @NotNull
        private String hgvsProtein = Strings.EMPTY;
        @NotNull
        private String consequence = Strings.EMPTY;
        @NotNull
        private String cosmicID = Strings.EMPTY;
        private int totalReadCount = 0;
        private int alleleReadCount = 0;

        public Builder() {
        }

        @NotNull
        public Builder gene(@NotNull final String gene) {
            this.gene = gene;
            return this;
        }

        @NotNull
        public Builder position(@NotNull final String position) {
            this.position = position;
            return this;
        }

        @NotNull
        public Builder ref(@NotNull final String ref) {
            this.ref = ref;
            return this;
        }

        @NotNull
        public Builder alt(@NotNull final String alt) {
            this.alt = alt;
            return this;
        }

        @NotNull
        public Builder transcript(@NotNull final String transcript) {
            this.transcript = transcript;
            return this;
        }

        @NotNull
        public Builder hgvsCoding(@NotNull final String hgvsCoding) {
            this.hgvsCoding = hgvsCoding;
            return this;
        }

        @NotNull
        public Builder hgvsProtein(@NotNull final String hgvsProtein) {
            this.hgvsProtein = hgvsProtein;
            return this;
        }

        @NotNull
        public Builder consequence(@NotNull final String consequence) {
            this.consequence = consequence;
            return this;
        }

        @NotNull
        public Builder cosmicID(@NotNull final String cosmicID) {
            this.cosmicID = cosmicID;
            return this;
        }

        @NotNull
        public Builder totalReadCount(final int totalReadCount) {
            this.totalReadCount = totalReadCount;
            return this;
        }

        @NotNull
        public Builder alleleReadCount(final int alleleReadCount) {
            this.alleleReadCount = alleleReadCount;
            return this;
        }

        @NotNull
        public VariantReport build() {
            return new VariantReport(gene, position, ref, alt, transcript, hgvsCoding, hgvsProtein, consequence,
                    cosmicID, totalReadCount, alleleReadCount);
        }
    }
}
