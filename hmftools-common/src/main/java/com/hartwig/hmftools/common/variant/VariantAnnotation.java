package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class VariantAnnotation {

    @NotNull
    private final String allele;
    @NotNull
    private final List<VariantConsequence> consequences;
    @NotNull
    private final String severity;
    @NotNull
    private final String gene;
    @NotNull
    private final String geneID;
    @NotNull
    private final String featureType;
    @NotNull
    private final String featureID;
    @NotNull
    private final String transcriptBioType;
    @NotNull
    private final String rank;
    @NotNull
    private final String hgvsCoding;
    @NotNull
    private final String hgvsProtein;
    @NotNull
    private final String cDNAPosAndLength;
    @NotNull
    private final String cdsPosAndLength;
    @NotNull
    private final String aaPosAndLength;
    @NotNull
    private final String distance;
    @NotNull
    private final String addition;

    private VariantAnnotation(@NotNull final String allele, @NotNull final List<VariantConsequence> consequences,
            @NotNull final String severity, @NotNull final String gene, @NotNull final String geneID,
            @NotNull final String featureType, @NotNull final String featureID,
            @NotNull final String transcriptBioType, @NotNull final String rank, @NotNull final String hgvsCoding,
            @NotNull final String hgvsProtein, @NotNull final String cDNAPosAndLength,
            @NotNull final String cdsPosAndLength, @NotNull final String aaPosAndLength,
            @NotNull final String distance, @NotNull final String addition) {
        this.allele = allele;
        this.consequences = consequences;
        this.severity = severity;
        this.gene = gene;
        this.geneID = geneID;
        this.featureType = featureType;
        this.featureID = featureID;
        this.transcriptBioType = transcriptBioType;
        this.rank = rank;
        this.hgvsCoding = hgvsCoding;
        this.hgvsProtein = hgvsProtein;
        this.cDNAPosAndLength = cDNAPosAndLength;
        this.cdsPosAndLength = cdsPosAndLength;
        this.aaPosAndLength = aaPosAndLength;
        this.distance = distance;
        this.addition = addition;
    }

    @NotNull
    String allele() {
        return allele;
    }

    @NotNull
    public List<VariantConsequence> consequences() {
        return consequences;
    }

    @NotNull
    String severity() {
        return severity;
    }

    @NotNull
    String gene() {
        return gene;
    }

    @NotNull
    String geneID() {
        return geneID;
    }

    @NotNull
    public String featureType() {
        return featureType;
    }

    @NotNull
    public String featureID() {
        return featureID;
    }

    @NotNull
    String transcriptBioType() {
        return transcriptBioType;
    }

    @NotNull
    String rank() {
        return rank;
    }

    @NotNull
    String hgvsCoding() {
        return hgvsCoding;
    }

    @NotNull
    String hgvsProtein() {
        return hgvsProtein;
    }

    @NotNull
    String cDNAPosAndLength() {
        return cDNAPosAndLength;
    }

    @NotNull
    String cdsPosAndLength() {
        return cdsPosAndLength;
    }

    @NotNull
    String aaPosAndLength() {
        return aaPosAndLength;
    }

    @NotNull
    String distance() {
        return distance;
    }

    @NotNull
    String addition() {
        return addition;
    }

    public static class Builder {
        @NotNull
        private String allele = Strings.EMPTY;
        @NotNull
        private List<VariantConsequence> consequences = Lists.newArrayList();
        @NotNull
        private String severity = Strings.EMPTY;
        @NotNull
        private String gene = Strings.EMPTY;
        @NotNull
        private String geneID = Strings.EMPTY;
        @NotNull
        private String featureType = Strings.EMPTY;
        @NotNull
        private String featureID = Strings.EMPTY;
        @NotNull
        private String transcriptBioType = Strings.EMPTY;
        @NotNull
        private String rank = Strings.EMPTY;
        @NotNull
        private String hgvsCoding = Strings.EMPTY;
        @NotNull
        private String hgvsProtein = Strings.EMPTY;
        @NotNull
        private String cDNAPosAndLength = Strings.EMPTY;
        @NotNull
        private String cdsPosAndLength = Strings.EMPTY;
        @NotNull
        private String aaPosAndLength = Strings.EMPTY;
        @NotNull
        private String distance = Strings.EMPTY;
        @NotNull
        private String addition = Strings.EMPTY;

        public Builder() {
        }

        @NotNull
        Builder allele(@NotNull final String allele) {
            this.allele = allele;
            return this;
        }

        @NotNull
        public Builder consequences(@NotNull final List<VariantConsequence> consequences) {
            this.consequences = consequences;
            return this;
        }

        @NotNull
        Builder severity(@NotNull final String severity) {
            this.severity = severity;
            return this;
        }

        @NotNull
        Builder gene(@NotNull final String gene) {
            this.gene = gene;
            return this;
        }

        @NotNull
        Builder geneID(@NotNull final String geneID) {
            this.geneID = geneID;
            return this;
        }

        @NotNull
        public Builder featureType(@NotNull final String featureType) {
            this.featureType = featureType;
            return this;
        }

        @NotNull
        public Builder featureID(@NotNull final String featureID) {
            this.featureID = featureID;
            return this;
        }

        @NotNull
        Builder transcriptBioType(@NotNull final String transcriptBioType) {
            this.transcriptBioType = transcriptBioType;
            return this;
        }

        @NotNull
        Builder rank(@NotNull final String rank) {
            this.rank = rank;
            return this;
        }

        @NotNull
        Builder hgvsCoding(@NotNull final String hgvsCoding) {
            this.hgvsCoding = hgvsCoding;
            return this;
        }

        @NotNull
        Builder hgvsProtein(@NotNull final String hgvsProtein) {
            this.hgvsProtein = hgvsProtein;
            return this;
        }

        @NotNull
        Builder cDNAPosAndLength(@NotNull final String cDNAPosAndLength) {
            this.cDNAPosAndLength = cDNAPosAndLength;
            return this;
        }

        @NotNull
        Builder cdsPosAndLength(@NotNull final String cdsPosAndLength) {
            this.cdsPosAndLength = cdsPosAndLength;
            return this;
        }

        @NotNull
        Builder aaPosAndLength(@NotNull final String aaPosAndLength) {
            this.aaPosAndLength = aaPosAndLength;
            return this;
        }

        @NotNull
        Builder distance(@NotNull final String distance) {
            this.distance = distance;
            return this;
        }

        @NotNull
        Builder addition(@NotNull final String addition) {
            this.addition = addition;
            return this;
        }

        @NotNull
        public VariantAnnotation build() {
            return new VariantAnnotation(allele, consequences, severity, gene, geneID, featureType, featureID,
                    transcriptBioType, rank, hgvsCoding, hgvsProtein, cDNAPosAndLength, cdsPosAndLength,
                    aaPosAndLength, distance, addition);
        }
    }
}
