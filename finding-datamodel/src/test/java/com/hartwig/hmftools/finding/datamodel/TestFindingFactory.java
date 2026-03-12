package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("unused")
public class TestFindingFactory
{

    @NotNull
    public static <T extends Finding> FindingList<T> buildFindingsList(@NotNull FindingsStatus findingsStatus, @NotNull List<T> findings)
    {
        return FindingListBuilder.<T>builder().status(findingsStatus).findings(findings).build();
    }

    @NotNull
    public static <T extends Driver> DriverFindingList<T> buildDriverFindingsList(@NotNull FindingsStatus findingsStatus,
            @NotNull List<T> findings)
    {
        return DriverFindingListBuilder.<T>builder().status(findingsStatus).findings(findings).build();
    }

    @NotNull
    public static <T> FindingItem<T> buildFindingItem(@NotNull FindingsStatus findingsStatus, @NotNull T finding)
    {
        return FindingItemBuilder.<T>builder().status(findingsStatus).finding(finding).build();
    }

    @NotNull
    public static PurityPloidyFitBuilder purityPloidyFitBuilder()
    {
        return PurityPloidyFitBuilder.builder()
                .fittedPurityMethod(PurityPloidyFit.FittedPurityMethod.NORMAL);
    }

    @NotNull
    public static QcBuilder qcBuilder()
    {
        return QcBuilder.builder();
    }

    @NotNull
    public static HomologousRecombinationBuilder homologousRecombinationBuilder()
    {
        return HomologousRecombinationBuilder.builder()
                .findingKey("")
                .hrdType("")
                .hrdValue(0)
                .hrStatus(HomologousRecombination.HrStatus.HR_DEFICIENT)
                .brca1Value(0)
                .brca2Value(0)
                .lohCopyNumbers(List.of())
                .relatedGenes(List.of());
    }

    @NotNull
    public static MicrosatelliteStabilityBuilder microsatelliteStabilityBuilder()
    {
        return MicrosatelliteStabilityBuilder.builder()
                .findingKey("")
                .microsatelliteStatus(MicrosatelliteStability.MicrosatelliteStatus.MSI)
                .lohCopyNumbers(List.of())
                .relatedGenes(List.of());
    }

    @NotNull
    public static TumorMutationalBurdenBuilder tumorMutationalBurdenBuilder()
    {
        return TumorMutationalBurdenBuilder.builder()
                .findingKey("")
                .status(TumorMutationalBurden.Status.HIGH)
                .burdenPerMb(0)
                .svBurden(0);
    }

    @NotNull
    public static TumorMutationalLoadBuilder tumorMutationalLoadBuilder()
    {
        return TumorMutationalLoadBuilder.builder()
                .findingKey("")
                .load(0)
                .status(TumorMutationalLoad.Status.HIGH);
    }

    @NotNull
    public static PredictedTumorOriginBuilder predictedTumorOriginBuilder()
    {
        return PredictedTumorOriginBuilder.builder()
                .findingKey("")
                .mode(PredictedTumorOrigin.CuppaMode.WGS)
                .cancerType("")
                .snvPairwiseClassifier(0.0)
                .genomicPositionClassifier(0.0)
                .featureClassifier(0.0)
                .altSjCohortClassifier(0.0)
                .expressionPairwiseClassifier(0.0);
    }

    @NotNull
    public static SmallVariantBuilder variantBuilder()
    {
        return SmallVariantBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .type(SmallVariant.VariantType.SNP)
                .gene("")
                .chromosome("")
                .position(0)
                .ref("")
                .alt("")
                .adjustedCopyNumber(0D)
                .minorAlleleCopyNumber(0D)
                .variantCopyNumber(0D)
                .hotspot(SmallVariant.HotspotType.NON_HOTSPOT)
                .allelicDepth(depthBuilder().build())
                .subclonalLikelihood(0D)
                .biallelic(false)
                .biallelicLikelihood(0D)
                .genotypeStatus(SmallVariant.GenotypeStatus.UNKNOWN)
                .worstCodingEffect(SmallVariant.CodingEffect.SPLICE)
                .adjustedVAF(0)
                .repeatCount(0)
                .transcriptImpact(transcriptImpactBuilder().reported(false).build());
    }

    @NotNull
    public static SmallVariantAllelicDepthBuilder depthBuilder()
    {
        return SmallVariantAllelicDepthBuilder.builder().totalReadCount(0).alleleReadCount(0);
    }

    @NotNull
    public static SmallVariantTranscriptImpactBuilder transcriptImpactBuilder()
    {
        return SmallVariantTranscriptImpactBuilder.builder()
                .transcript("")
                .hgvsCodingImpact("")
                .hgvsProteinImpact("")
                .inSpliceRegion(false)
                .codingEffect(SmallVariant.CodingEffect.UNDEFINED)
                .effects(Set.of())
                .reported(true);
    }

    @NotNull
    public static DisruptionBuilder disruptionBuilder()
    {
        return DisruptionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .gene("")
                .chromosome("")
                .chromosomeBand("")
                .breakendType(Disruption.BreakendType.BND)
                .transcript("")
                .isCanonical(true)
                .breakendStart(breakendBuilder().build())
                .breakendEnd(breakendBuilder().build())
                .type(Disruption.Type.DISRUPTION);
    }

    @NotNull
    public static BreakendBuilder breakendBuilder()
    {
        return BreakendBuilder.builder()
                .geneOrientation(Breakend.GeneOrientation.UPSTREAM)
                .regionType(Breakend.TranscriptRegionType.DOWNSTREAM)
                .codingType(Breakend.TranscriptCodingType.CODING);
    }

    @NotNull
    public static DisruptionBuilder homozygousDisruptionBuilder()
    {
        return DisruptionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .chromosomeBand("")
                .chromosome("")
                .transcript("")
                .isCanonical(true)
                .breakendType(Disruption.BreakendType.DEL)
                .type(Disruption.Type.HOM_DEL_DISRUPTION);
    }

    @NotNull
    public static GainDeletionBuilder gainDeletionBuilder()
    {
        return GainDeletionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .chromosome("1")
                .chromosomeBand("p")
                .gene("")
                .transcript("")
                .isCanonical(false)
                .somaticType(GainDeletion.Type.HOM_DEL)
                .geneExtent(GainDeletion.GeneExtent.FULL_GENE)
                .tumorMinCopyNumber(0)
                .tumorMaxCopyNumber(0)
                .tumorMinMinorAlleleCopyNumber(0)
                .chromosomeArmCopyNumber(0);
    }

    @NotNull
    public static FusionBuilder fusionBuilder()
    {
        return FusionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .reportedType(Fusion.FusionType.NONE)
                .geneStart("")
                .geneTranscriptStart("")
                .geneContextStart("")
                .fusedExonUp(0)
                .geneEnd("")
                .geneTranscriptEnd("")
                .geneContextEnd("")
                .fusedExonDown(0)
                .phased(Fusion.FusionPhasedType.OUT_OF_FRAME)
                .junctionCopyNumber(0D)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept(List.of())
                .domainsLost(List.of());
    }

    public static <T extends Driver> DriverFindingListBuilder<T> driverFindingsBuilder(List<T> findings)
    {
        return DriverFindingListBuilder.<T>builder().status(FindingsStatus.OK).findings(findings);
    }

    @NotNull
    public static VirusBuilder virusBuilder(boolean reported, DriverInterpretation driverInterpretation)
    {
        return VirusBuilder.builder()
                .driver(driverFields(reported, driverInterpretation))
                .name("")
                .qcStatus(Virus.VirusBreakendQCStatus.NO_ABNORMALITIES);
    }

    @NotNull
    public static HlaAlleleBuilder hlaAlleleBuilder()
    {
        return HlaAlleleBuilder.builder()
                .findingKey("")
                .gene("")
                .geneClass("")
                .germlineCopyNumber(0)
                .alleleGroup("")
                .hlaProtein("")
                .tumorCopyNumber(0D)
                .refFragments(0)
                .tumorFragments(0)
                .rnaFragments(0)
                .somaticMissense(0D)
                .somaticNonsenseOrFrameshift(0D)
                .somaticSplice(0D)
                .somaticSynonymous(0D)
                .somaticInframeIndel(0D);
    }

    public static DriverFields driverFields(boolean reported, DriverInterpretation interpretation)
    {
        return driverFieldsBuilder()
                .driverInterpretation(interpretation)
                .reportedStatus(reported ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED)
                .build();
    }

    public static DriverFieldsBuilder driverFieldsBuilder()
    {
        return DriverFieldsBuilder.builder()
                .findingKey("")
                .driverSource(DriverSource.SOMATIC);
    }
}
