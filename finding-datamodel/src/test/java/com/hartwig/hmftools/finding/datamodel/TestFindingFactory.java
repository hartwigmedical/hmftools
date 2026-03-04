package com.hartwig.hmftools.finding.datamodel;

import java.util.List;
import java.util.Objects;
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
                .qc(TestFindingFactory.qcBuilder().status(Set.of()).build())
                .fittedPurityMethod(PurityPloidyFit.FittedPurityMethod.NORMAL);
    }

    @NotNull
    public static PurityPloidyFitQcBuilder qcBuilder()
    {
        return PurityPloidyFitQcBuilder.builder();
    }

    @NotNull
    public static HomologousRecombinationBuilder homologousRecombinationBuilder()
    {
        return HomologousRecombinationBuilder.builder()
                .findingKey("")
                .hrdType("")
                .hrdValue(0)
                .hrStatus(HomologousRecombination.HrStatus.UNKNOWN)
                .brca1Value(0)
                .brca2Value(0)
                .lohCopyNumbers(List.of())
                .genes(List.of());
    }

    @NotNull
    public static MicrosatelliteStabilityBuilder microsatelliteStabilityBuilder()
    {
        return MicrosatelliteStabilityBuilder.builder()
                .findingKey("")
                .microsatelliteStatus(MicrosatelliteStability.MicrosatelliteStatus.UNKNOWN)
                .lohCopyNumbers(List.of())
                .genes(List.of());
    }

    @NotNull
    public static TumorMutationStatusBuilder mutationStatusBuilder()
    {
        return TumorMutationStatusBuilder.builder()
                .findingKey("")
                .tumorMutationalLoad(0)
                .tumorMutationalBurdenStatus(TumorMutationStatus.Status.UNKNOWN)
                .tumorMutationalBurdenPerMb(0)
                .svTumorMutationalBurden(0)
                .tumorMutationalLoadStatus(TumorMutationStatus.Status.UNKNOWN);
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
                .biallelicProbability(0D)
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
                .breakendType(Breakend.Type.BND)
                .transcript("")
                .isCanonical(true)
                .breakendStart(breakendBuilder().build())
                .breakendEnd(breakendBuilder().build())
                .type(Disruption.Type.SOMATIC_DISRUPTION);
    }

    @NotNull
    public static BreakendBuilder breakendBuilder()
    {
        return BreakendBuilder.builder()
                .gene("")
                .chromosome("")
                .chromosomeBand("")
                .transcript("")
                .geneOrientation(Breakend.GeneOrientation.UPSTREAM)
                .type(Breakend.Type.BND)
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
                .breakendType(Breakend.Type.DEL)
                .type(Disruption.Type.SOMATIC_HOM_DEL_DISRUPTION);
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
                .tumorMinCopies(0)
                .tumorMaxCopies(0)
                .tumorMinMinorAlleleCopies(0)
                .chromosomeArmCopies(0);
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
                .domainsKept("")
                .domainsLost("");
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
                .event("")
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
        return DriverFieldsBuilder.builder()
                .findingKey("")
                .event("")
                .driverSource(DriverSource.SOMATIC)
                .driverInterpretation(interpretation)
                .reportedStatus(reported ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED)
                .build();
    }
}
