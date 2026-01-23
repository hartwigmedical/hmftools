package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

@SuppressWarnings("unused")
public class TestFindingFactory {

    @NotNull
    public static <T extends Finding> FindingList<T> buildFindingsList(@NotNull FindingsStatus findingsStatus, @NotNull List<T> findings) {
        return FindingListBuilder.<T>builder().status(findingsStatus).findings(findings).build();
    }

    @NotNull
    public static <T extends Driver> DriverFindingList<T> buildDriverFindingsList(@NotNull FindingsStatus findingsStatus,
            @NotNull List<T> findings) {
        return DriverFindingListBuilder.<T>builder().status(findingsStatus).findings(findings).build();
    }

    @NotNull
    public static <T> FindingItem<T> buildFindingItem(@NotNull FindingsStatus findingsStatus, @NotNull T finding) {
        return FindingItemBuilder.<T>builder().status(findingsStatus).finding(finding).build();
    }

    @NotNull
    public static PurityPloidyFitBuilder purityPloidyFitBuilder() {
        return PurityPloidyFitBuilder.builder().qc(TestFindingFactory.qcBuilder().status(Set.of()).build());
    }

    @NotNull
    public static PurityPloidyFitQcBuilder qcBuilder() {
        return PurityPloidyFitQcBuilder.builder();
    }

    @NotNull
    public static HomologousRecombinationBuilder homologousRecombinationBuilder() {
        return HomologousRecombinationBuilder.builder()
                .findingKey("")
                .hrdType("")
                .hrdValue(0)
                .hrStatus(ChordStatus.UNKNOWN)
                .brca1Value(0)
                .brca2Value(0)
                .lohCopyNumbers(List.of())
                .genes(List.of());
    }

    @NotNull
    public static MicrosatelliteStabilityBuilder microsatelliteStabilityBuilder() {
        return MicrosatelliteStabilityBuilder.builder()
                .findingKey("")
                .microsatelliteStatus(PurpleMicrosatelliteStatus.UNKNOWN)
                .lohCopyNumbers(List.of())
                .genes(List.of());
    }

    @NotNull
    public static TumorMutationStatusBuilder mutationStatusBuilder() {
        return TumorMutationStatusBuilder.builder()
                .findingKey("")
                .tumorMutationalLoad(0)
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(0)
                .svTumorMutationalBurden(0)
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.UNKNOWN);
    }

    @NotNull
    public static PredictedTumorOriginBuilder predictedTumorOriginBuilder() {
        return PredictedTumorOriginBuilder.builder().findingKey("");
    }

    @NotNull
    public static SmallVariantBuilder variantBuilder() {
        return SmallVariantBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .type(PurpleVariantType.SNP)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .adjustedCopyNumber(0D)
                .minorAlleleCopyNumber(0D)
                .variantCopyNumber(0D)
                .hotspot(HotspotType.NON_HOTSPOT)
                .allelicDepth(depthBuilder().build())
                .subclonalLikelihood(0D)
                .biallelic(false)
                .biallelicProbability(0D)
                .genotypeStatus(PurpleGenotypeStatus.UNKNOWN)
                .worstCodingEffect(PurpleCodingEffect.SPLICE)
                .adjustedVAF(0)
                .repeatCount(0)
                .transcriptImpact(transcriptImpactBuilder().reported(false).build());
    }

    @NotNull
    public static SmallVariantAllelicDepthBuilder depthBuilder() {
        return SmallVariantAllelicDepthBuilder.builder().totalReadCount(0).alleleReadCount(0);
    }

    @NotNull
    public static SmallVariantTranscriptImpactBuilder transcriptImpactBuilder() {
        return SmallVariantTranscriptImpactBuilder.builder()
                .transcript(Strings.EMPTY)
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .inSpliceRegion(false)
                .codingEffect(PurpleCodingEffect.UNDEFINED)
                .effects(Set.of())
                .reported(true);
    }

    @NotNull
    public static DisruptionBuilder disruptionBuilder() {
        return DisruptionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .breakendType(LinxBreakendType.BND)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .breakendStart(breakendBuilder().build())
                .breakendEnd(breakendBuilder().build())
                .type(Disruption.Type.SOMATIC_DISRUPTION);
    }

    @NotNull
    public static BreakendBuilder breakendBuilder() {
        return BreakendBuilder.builder();
    }

    @NotNull
    public static DisruptionBuilder homozygousDisruptionBuilder() {
        return DisruptionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .chromosomeBand("")
                .chromosome("")
                .transcript("")
                .isCanonical(true)
                .breakendType(LinxBreakendType.DEL)
                .type(Disruption.Type.SOMATIC_HOM_DEL_DISRUPTION);
    }

    @NotNull
    public static GainDeletionBuilder gainDeletionBuilder() {
        return GainDeletionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .chromosome("1")
                .chromosomeBand("p")
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false)
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .tumorMinCopies(0)
                .tumorMaxCopies(0)
                .tumorMinMinorAlleleCopies(0)
                .chromosomeArmCopies(0)
                .type(GainDeletion.Type.SOMATIC_GAIN);
    }

    @NotNull
    public static FusionBuilder fusionBuilder() {
        return FusionBuilder.builder()
                .driver(driverFields(true, DriverInterpretation.HIGH))
                .reportedType(LinxFusionType.NONE)
                .geneStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .fusedExonUp(0)
                .geneEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .fusedExonDown(0)
                .phased(FusionPhasedType.OUT_OF_FRAME)
                .junctionCopyNumber(0D)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY);
    }

    public static <T extends Driver> DriverFindingListBuilder<T> driverFindingsBuilder(List<T> findings) {
        return DriverFindingListBuilder.<T>builder().status(FindingsStatus.OK).findings(findings);
    }

    @NotNull
    public static VirusBuilder virusBuilder(boolean reported, DriverInterpretation driverInterpretation) {
        return VirusBuilder.builder().driver(driverFields(reported, driverInterpretation));
    }

    @NotNull
    public static HlaAlleleBuilder hlaAlleleBuilder() {
        return HlaAlleleBuilder.builder()
                .findingKey("")
                .gene("")
                .germlineCopyNumber(0)
                .allele(Strings.EMPTY)
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

    public static DriverFields driverFields(boolean reported, DriverInterpretation interpretation) {
        return DriverFieldsBuilder.builder()
                .findingKey("")
                .driverSource(DriverSource.SOMATIC)
                .driverInterpretation(interpretation)
                .reportedStatus(reported ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED)
                .build();
    }
}
