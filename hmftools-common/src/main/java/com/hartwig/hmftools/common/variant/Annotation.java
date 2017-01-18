package com.hartwig.hmftools.common.variant;

import org.jetbrains.annotations.NotNull;

class Annotation {

    @NotNull
    private final String allele;
    @NotNull
    private final String annotation;
    @NotNull
    private final String annotationImpact;
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

    Annotation(@NotNull final String allele, @NotNull final String annotation, @NotNull final String annotationImpact,
            @NotNull final String gene, @NotNull final String geneID, @NotNull final String featureType,
            @NotNull final String featureID, @NotNull final String transcriptBioType, @NotNull final String rank,
            @NotNull final String hgvsCoding, @NotNull final String hgvsProtein,
            @NotNull final String cDNAPosAndLength, @NotNull final String cdsPosAndLength,
            @NotNull final String aaPosAndLength, @NotNull final String distance, @NotNull final String addition) {
        this.allele = allele;
        this.annotation = annotation;
        this.annotationImpact = annotationImpact;
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
    String annotation() {
        return annotation;
    }

    @NotNull
    String annotationImpact() {
        return annotationImpact;
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
    String featureType() {
        return featureType;
    }

    @NotNull
    String featureID() {
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
}
