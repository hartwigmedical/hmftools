package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BRCApart2 {

    @NotNull
    public abstract String alleleCountFINExAC();

    @NotNull
    public abstract String conditionCategoryENIGMA();

    @NotNull
    public abstract String alleleFrequencyESP();

    @NotNull
    public abstract String homozygousCountOTHExAC();

    @NotNull
    public abstract String geneticOriginLOVD();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String homozygousCountAMRExAC();

    @NotNull
    public abstract String clinicalSignificanceClinVar();

    @NotNull
    public abstract String aaAlleleFrequencyESP();

    @NotNull
    public abstract String proteinChange();

    @NotNull
    public abstract String variantInExLOVD();

    @NotNull
    public abstract String eaAlleleFrequencyESP();

    @NotNull
    public abstract String hgvsRNA();

    @NotNull
    public abstract String clinicalSignificanceCitationsENIGMA();

    @NotNull
    public abstract String variantEffectLOVD();

    @NotNull
    public abstract String polyphenPrediction();

    @NotNull
    public abstract String dataReleaseId();

    @NotNull
    public abstract String hg37Start();

    @NotNull
    public abstract String hg36Start();

    @NotNull
    public abstract String siftScore();

    @NotNull
    public abstract String genomicCoordinateHg38();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String literatureCitationBIC();

    @NotNull
    public abstract String variantHaplotypeLOVD();

    @NotNull
    public abstract String alleleFrequencyNFEExAC();

    @NotNull
    public abstract String hg38Start();

    @NotNull
    public abstract String pos();

    @NotNull
    public abstract String dateLastEvaluatedENIGMA();

    @NotNull
    public abstract String alleleNumberSASExAC();

    @NotNull
    public abstract String alleleNumberAMRExAC();

    @NotNull
    public abstract String dbIdLOVD();
}
