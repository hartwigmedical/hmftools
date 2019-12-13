package com.hartwig.hmftools.vicc.datamodel.brca;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BRCA implements KbSpecificObject {

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String pos();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String genomicCoordinateHg36();

    @NotNull
    public abstract String hg36Start();

    @NotNull
    public abstract String hg36End();

    @NotNull
    public abstract String genomicCoordinateHg37();

    @NotNull
    public abstract String hg37Start();

    @NotNull
    public abstract String hg37End();

    @NotNull
    public abstract String genomicCoordinateHg38();

    @NotNull
    public abstract String hg38Start();

    @NotNull
    public abstract String hg38End();

    @NotNull
    public abstract String proteinChange();

    @NotNull
    public abstract String referenceSequence();

    @NotNull
    public abstract String synonyms();

    @NotNull
    public abstract String hgvsCDNA();

    @NotNull
    public abstract String hgvsProtein();

    @NotNull
    public abstract String hgvsRNA();

    @NotNull
    public abstract String siftScore();

    @NotNull
    public abstract String siftPrediction();

    @NotNull
    public abstract String polyphenScore();

    @NotNull
    public abstract String polyphenPrediction();

    @NotNull
    public abstract String pathogenicityAll();

    @NotNull
    public abstract String pathogenicityExpert();

    @NotNull
    public abstract String alleleFrequency();

    @NotNull
    public abstract String maxAlleleFrequency();

    @NotNull
    public abstract String discordant();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String changeTypeId();

    @NotNull
    public abstract String dataReleaseId();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String sourceURL();

    @NotNull
    public abstract BRCAAnnotation1000Genomes annotation1000Genomes();

    @NotNull
    public abstract BRCAAnnotationBIC annotationBIC();

    @NotNull
    public abstract BRCAAnnotationClinVar annotationClinVar();

    @NotNull
    public abstract BRCAAnnotationENIGMA annotationENIGMA();

    @NotNull
    public abstract BRCAAnnotationESP annotationESP();

    @NotNull
    public abstract BRCAAnnotationExAC annotationExAC();

    @NotNull
    public abstract BRCAAnnotationExLOVD annotationExLOVD();

    @NotNull
    public abstract BRCAAnnotationLOVD annotationLOVD();

}
