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
    @NotNull
    private final String alleleFrequency;
    @NotNull
    private final String readDepth;

    private VariantReport(@NotNull final String gene, @NotNull final String position, @NotNull final String ref,
            @NotNull final String alt, @NotNull final String transcript, @NotNull final String hgvsCoding,
            @NotNull final String hgvsProtein, @NotNull final String consequence, @NotNull final String cosmicID,
            @NotNull final String alleleFrequency, @NotNull final String readDepth) {
        this.gene = gene;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.transcript = transcript;
        this.hgvsCoding = hgvsCoding;
        this.hgvsProtein = hgvsProtein;
        this.consequence = consequence;
        this.cosmicID = cosmicID;
        this.alleleFrequency = alleleFrequency;
        this.readDepth = readDepth;
    }

    @NotNull
    public String getGene() {
        return gene;
    }

    @NotNull
    public String getPosition() {
        return position;
    }

    @NotNull
    public String getRef() {
        return ref;
    }

    @NotNull
    public String getAlt() {
        return alt;
    }

    @NotNull
    public String getTranscript() {
        return transcript;
    }

    @NotNull
    public String getHgvsCoding() {
        return hgvsCoding;
    }

    @NotNull
    public String getHgvsProtein() {
        return hgvsProtein;
    }

    @NotNull
    public String getConsequence() {
        return consequence;
    }

    @NotNull
    public String getCosmicID() {
        return cosmicID;
    }

    @NotNull
    public String getAlleleFrequency() {
        return alleleFrequency;
    }

    @NotNull
    public String getReadDepth() {
        return readDepth;
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
        @NotNull
        private String alleleFrequency = Strings.EMPTY;
        @NotNull
        private String readDepth = Strings.EMPTY;

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
        public Builder alleleFrequency(@NotNull final String alleleFrequency) {
            this.alleleFrequency = alleleFrequency;
            return this;
        }

        @NotNull
        public Builder readDepth(@NotNull final String readDepth) {
            this.readDepth = readDepth;
            return this;
        }

        @NotNull
        public VariantReport build() {
            return new VariantReport(gene, position, ref, alt, transcript, hgvsCoding, hgvsProtein, consequence,
                    cosmicID, alleleFrequency, readDepth);
        }
    }
}
