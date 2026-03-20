package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.DisruptionFactory.createGermlineDisruptions;
import static com.hartwig.hmftools.finding.DisruptionFactory.createSomaticDisruptions;
import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.ResultIssue.REF_REQUIRED;
import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.ResultIssue.WGS_REQUIRED;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangePlots;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.Genes;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.FusionBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombination;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombinationBuilder;
import com.hartwig.hmftools.finding.datamodel.MetaPropertiesBuilder;
import com.hartwig.hmftools.finding.datamodel.MicrosatelliteStability;
import com.hartwig.hmftools.finding.datamodel.MicrosatelliteStabilityBuilder;
import com.hartwig.hmftools.finding.datamodel.PharmacoGenotype;
import com.hartwig.hmftools.finding.datamodel.PharmacoGenotypeBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOrigin;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginPredictionBuilder;
import com.hartwig.hmftools.finding.datamodel.PurityPloidyFit;
import com.hartwig.hmftools.finding.datamodel.PurityPloidyFitBuilder;
import com.hartwig.hmftools.finding.datamodel.Qc;
import com.hartwig.hmftools.finding.datamodel.QcBuilder;
import com.hartwig.hmftools.finding.datamodel.RefGenomeVersion;
import com.hartwig.hmftools.finding.datamodel.SequencingScope;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.TumorMutationalBurden;
import com.hartwig.hmftools.finding.datamodel.TumorMutationalBurdenBuilder;
import com.hartwig.hmftools.finding.datamodel.TumorMutationalLoad;
import com.hartwig.hmftools.finding.datamodel.TumorMutationalLoadBuilder;
import com.hartwig.hmftools.finding.datamodel.Virus;
import com.hartwig.hmftools.finding.datamodel.VirusBuilder;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory
{
    public static FindingRecord fromOrangeJsonWithTranscriptFile(Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, Gender gender) throws IOException
    {
        try(Reader reader = Files.newBufferedReader(orangeJson))
        {
            OrangeRecord orangeRecord = OrangeJson.getInstance().read(reader);
            return fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv, gender);
        }
    }

    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, Gender gender) throws IOException
    {
        FindingConfig findingConfig =
                FindingConfig.createFindingConfig(clinicalTranscriptsTsv, driverGeneTsv, orangeRecord.refGenomeVersion(), gender);

        LinxRecord linx = orangeRecord.linx();
        PurpleRecord purple = orangeRecord.purple();

        FindingStatus findingStatus = FindingsStatusFactory.toFindingsStatus(purple.fit().qc().status());

        boolean hasRefSample = orangeRecord.refSample() != null;

        ExperimentType experimentType = orangeRecord.experimentType();
        FindingRecordBuilder
                builder = FindingRecordBuilder.builder()
                .version("1.0")
                .metaProperties(MetaPropertiesBuilder.builder()
                        .refGenomeVersion(RefGenomeVersion.valueOf(orangeRecord.refGenomeVersion().name()))
                        .sequencingScope(SequencingScope.valueOf(experimentType.name()))
                        .pipelineVersion(orangeRecord.pipelineVersion())
                        .sampleId(orangeRecord.sampleId())
                        .samplingDate(orangeRecord.samplingDate())
                        .build())
                .qc(createQc(purple))
                .fusions(createFusionsFindings(orangeRecord.linx(), findingStatus));

        DriverFindingList<GainDeletion>
                somaticGainDeletions = addPurpleFindings(builder, orangeRecord, findingStatus, findingConfig);

        builder.somaticDisruptions(createSomaticDisruptions(linx, findingStatus))
                .germlineDisruptions(createGermlineDisruptions(orangeRecord.refSample() != null, linx, findingStatus))
                .viruses(createVirusFindings(orangeRecord.virusInterpreter(), experimentType, findingStatus))
                .homologousRecombination(createHomologousRecombination(orangeRecord.chord(), purple, linx, somaticGainDeletions, findingStatus, experimentType, hasRefSample));

        return builder.predictedTumorOrigin(createPredictedTumorOrigin(orangeRecord.cuppa(), orangeRecord.plots(), experimentType, findingStatus))
                .hlaAlleles(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, findingStatus))
                .pharmacoGenotypes(createPharmacoGenotypesFindings(orangeRecord.peach(), findingStatus))
                .build();
    }

    // return the gain deletions cause they are needed by HRD, will see if we can find a better way
    private static DriverFindingList<GainDeletion> addPurpleFindings(
            FindingRecordBuilder builder, final OrangeRecord orangeRecord,
            FindingStatus findingStatus,
            FindingConfig findingConfig)
    {
        boolean hasRefSample = orangeRecord.refSample() != null;

        PurpleRecord purple = orangeRecord.purple();

        builder.purityPloidyFit(createPurityPloidyFit(purple, orangeRecord.plots()));

        DriverFindingList<GainDeletion> somaticGainDeletions;

        DriverFindingList<SmallVariant> smallVariants =
                SmallVariantFactory.somaticSmallVariantFindings(purple, findingStatus, findingConfig);
        somaticGainDeletions =
                GainDeletionFactory.somaticGainDeletionFindings(orangeRecord.refGenomeVersion(), findingStatus, purple, findingConfig.gender());

        builder.somaticSmallVariants(smallVariants)
                .germlineSmallVariants(SmallVariantFactory.germlineSmallVariantFindings(hasRefSample, purple, findingStatus, findingConfig))
                .somaticGainDeletions(somaticGainDeletions)
                .germlineGainDeletions(GainDeletionFactory.germlineGainDeletionFindings(hasRefSample, orangeRecord.refGenomeVersion(), findingStatus, purple, findingConfig.gender()))
                .microsatelliteStability(createMicrosatelliteStability(purple, orangeRecord.linx(), somaticGainDeletions, findingStatus))
                .tumorMutationalLoad(createTumorMutationalLoad(purple, findingStatus))
                .tumorMutationalBurden(createTumorMutationalBurden(purple, findingStatus))
                .chromosomeArmCopyNumbers(new ArmCopyNumberFactory(purple.allSomaticCopyNumbers(), purple.fit().ploidy(),
                        findingConfig.gender(), orangeRecord.refGenomeVersion()).toArmCopyNumberFindings(findingStatus));

        return somaticGainDeletions;
    }

    private static Qc createQc(PurpleRecord purple)
    {
        PurpleFit purpleFit = purple.fit();
        Set<Qc.QCStatus> qcStatuses = purpleFit.qc().status().stream()
                .map(o -> switch(o)
                {
                    case PASS -> Qc.QCStatus.PASS;
                    case WARN_DELETED_GENES -> Qc.QCStatus.WARN_DELETED_GENES;
                    case WARN_HIGH_COPY_NUMBER_NOISE -> Qc.QCStatus.WARN_HIGH_COPY_NUMBER_NOISE;
                    case WARN_GENDER_MISMATCH -> Qc.QCStatus.WARN_GENDER_MISMATCH;
                    case WARN_LOW_PURITY -> Qc.QCStatus.WARN_LOW_PURITY;
                    case WARN_TINC -> Qc.QCStatus.WARN_TUMOR_IN_NORMAL_CONTAMINATION;
                    case FAIL_CONTAMINATION -> Qc.QCStatus.FAIL_CONTAMINATION;
                    case FAIL_NO_TUMOR -> Qc.QCStatus.FAIL_NO_TUMOR;
                    case FAIL_TINC -> Qc.QCStatus.FAIL_TUMOR_IN_NORMAL_CONTAMINATION;
                })
                .collect(Collectors.toSet());
        Set<Qc.GermlineAberration> germlineAberrations = purpleFit.qc().germlineAberrations().stream()
                .map(o -> switch(o)
                {
                    case NONE -> Qc.GermlineAberration.NONE;
                    case MOSAIC_X -> Qc.GermlineAberration.MOSAIC_X;
                    case KLINEFELTER -> Qc.GermlineAberration.KLINEFELTER;
                    case XYY -> Qc.GermlineAberration.XYY;
                    case TRISOMY_X -> Qc.GermlineAberration.TRISOMY_X;
                    case TRISOMY_13 -> Qc.GermlineAberration.TRISOMY_13;
                    case TRISOMY_15 -> Qc.GermlineAberration.TRISOMY_15;
                    case TRISOMY_18 -> Qc.GermlineAberration.TRISOMY_18;
                    case TRISOMY_21 -> Qc.GermlineAberration.TRISOMY_21;
                })
                .collect(Collectors.toSet());

        return QcBuilder.builder()
                .status(qcStatuses)
                .germlineAberrations(germlineAberrations)
                .amberMeanDepth(purpleFit.qc().amberMeanDepth())
                .contamination(purpleFit.qc().contamination())
                .totalCopyNumberSegments(purpleFit.qc().totalCopyNumberSegments())
                .unsupportedCopyNumberSegments(purpleFit.qc().unsupportedCopyNumberSegments())
                .deletedGenes(purpleFit.qc().deletedGenes())
                .build();
    }

    private static PurityPloidyFit createPurityPloidyFit(PurpleRecord purple, OrangePlots orangePlots)
    {
        PurpleFit purpleFit = purple.fit();

        return PurityPloidyFitBuilder.builder()
                .fittedPurityMethod(PurityPloidyFit.FittedPurityMethod.valueOf(purpleFit.fittedPurityMethod().name()))
                .purity(purpleFit.purity())
                .minPurity(purpleFit.minPurity())
                .maxPurity(purpleFit.maxPurity())
                .ploidy(purpleFit.ploidy())
                .minPloidy(purpleFit.minPloidy())
                .maxPloidy(purpleFit.maxPloidy())
                .purpleInputPlot(VisualisationFileUtil.create(orangePlots.purpleInputPlot()))
                .purpleCircosPlot(VisualisationFileUtil.create(orangePlots.purpleFinalCircosPlot()))
                .purpleClonalityPlot(VisualisationFileUtil.create(orangePlots.purpleClonalityPlot()))
                .purpleCopyNumberPlot(VisualisationFileUtil.create(orangePlots.purpleCopyNumberPlot()))
                .purpleVariantCopyNumberPlot(VisualisationFileUtil.create(orangePlots.purpleVariantCopyNumberPlot()))
                .purplePurityRangePlot(VisualisationFileUtil.create(orangePlots.purplePurityRangePlot()))
                .purpleKataegisPlot(VisualisationFileUtil.create(orangePlots.purpleKataegisPlot()))
                .build();
    }

    private static FindingItem<PredictedTumorOrigin> createPredictedTumorOrigin(@Nullable CuppaData cuppa, OrangePlots orangePlots,
            @NotNull ExperimentType experimentType,
            @NotNull FindingStatus findingStatus)
    {
        if(cuppa != null)
        {
            return FindingItemBuilder.<PredictedTumorOrigin>builder()
                    .status(findingStatus)
                    .finding(PredictedTumorOriginBuilder.builder()
                            .findingKey("predictedTumorOrigin")
                            .mode(cuppaMode(cuppa.mode()))
                            .predictions(cuppa.predictions().stream()
                                    .map(FindingRecordFactory::createPredictedTumorOriginPrediction)
                                    .toList())
                            .visualisationFile(VisualisationFileUtil.createNullable(orangePlots.cuppaSummaryPlot()))
                            .build()
                    )
                    .build();
        }
        else
        {
            if(experimentType == ExperimentType.TARGETED)
            {
                return FindingUtil.notAvailableFindingItem(Set.of(WGS_REQUIRED));
            }
            else
            {
                throw new IllegalStateException("Cuppa should not be null");
            }
        }
    }

    private static PredictedTumorOrigin.Prediction createPredictedTumorOriginPrediction(CuppaPrediction prediction)
    {
        return PredictedTumorOriginPredictionBuilder.builder()
                .findingKey(FindingKeys.predictedTumorOrigin(prediction.cancerType()))
                .cancerType(prediction.cancerType())
                .likelihood(prediction.likelihood())
                .snvPairwiseClassifier(prediction.snvPairwiseClassifier())
                .genomicPositionClassifier(prediction.genomicPositionClassifier())
                .featureClassifier(prediction.featureClassifier())
                .altSjCohortClassifier(prediction.altSjCohortClassifier())
                .expressionPairwiseClassifier(prediction.expressionPairwiseClassifier())
                .build();
    }

    private static PredictedTumorOrigin.CuppaMode cuppaMode(CuppaMode cuppaMode)
    {
        return switch(cuppaMode)
        {
            case WGS -> PredictedTumorOrigin.CuppaMode.WGS;
            case WGTS -> PredictedTumorOrigin.CuppaMode.WGTS;
        };
    }

    private static FindingItem<TumorMutationalBurden> createTumorMutationalBurden(PurpleRecord purple, FindingStatus findingStatus)
    {
        TumorMutationalBurden.Status status = tumorMutationalBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus());
        return FindingItemBuilder.<TumorMutationalBurden>builder()
                .status(findingStatus)
                .finding(TumorMutationalBurdenBuilder.builder()
                        .findingKey(FindingKeys.tumorMutationBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus()))
                        .burdenPerMb(ThresholdValueFactory.tmbValue(purple.characteristics().tumorMutationalBurdenPerMb()))
                        .status(status)
                        .svBurden(purple.characteristics().svTumorMutationalBurden())
                        .build())
                .build();
    }

    @NotNull
    private static TumorMutationalBurden.Status tumorMutationalBurdenStatus(PurpleTumorMutationalStatus status)
    {
        return switch(status)
        {
            case HIGH -> TumorMutationalBurden.Status.HIGH;
            case LOW -> TumorMutationalBurden.Status.LOW;
            case UNKNOWN -> throw new IllegalStateException("Tumor mutational burden status should not be UNKNOWN");
        };
    }

    private static FindingItem<TumorMutationalLoad> createTumorMutationalLoad(PurpleRecord purple, FindingStatus findingStatus)
    {
        TumorMutationalLoad.Status status = tumorMutationalLoadStatus(purple.characteristics().tumorMutationalLoadStatus());
        return FindingItemBuilder.<TumorMutationalLoad>builder()
                .status(findingStatus)
                .finding(TumorMutationalLoadBuilder.builder()
                        .findingKey(FindingKeys.tumorMutationLoadStatus(purple.characteristics().tumorMutationalLoadStatus()))
                        .load(ThresholdValueFactory.tmlValue(purple.characteristics().tumorMutationalLoad()))
                        .status(status)
                        .build())
                .build();
    }

    @NotNull
    private static TumorMutationalLoad.Status tumorMutationalLoadStatus(PurpleTumorMutationalStatus status)
    {
        return switch(status)
        {
            case HIGH -> TumorMutationalLoad.Status.HIGH;
            case LOW -> TumorMutationalLoad.Status.LOW;
            case UNKNOWN -> throw new IllegalStateException("Tumor mutational burden load should not be UNKNOWN");
        };
    }

    private static FindingItem<HomologousRecombination> createHomologousRecombination(@Nullable ChordRecord chord,
            PurpleRecord purple,
            LinxRecord linx,
            DriverFindingList<GainDeletion> gainDeletions,
            FindingStatus findingStatus,
            @NotNull ExperimentType experimentType,
            boolean hasRefSample)
    {
        if(chord != null)
        {
            HomologousRecombination.Status hrStatus = hrStatus(chord);
            boolean isPresent = hrStatus == HomologousRecombination.Status.HR_DEFICIENT;
            List<GainDeletion> lohGainDeletions = isPresent
                    ? filterLohGainDeletions(gainDeletions, Genes.HRD_GENES)
                    : List.of();
            List<String> drivingGenes = isPresent ? GeneListUtil.genes(purple.reportableSomaticVariants(),
                    purple.reportableSomaticGainsDels(),
                    linx.germlineHomozygousDisruptions(),
                    Genes.HRD_GENES) : List.of();
            return FindingItemBuilder.<HomologousRecombination>builder()
                    .status(findingStatus)
                    .finding(HomologousRecombinationBuilder.builder()
                            .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                            .status(hrStatus)
                            .hrdValue(ThresholdValueFactory.hrdValue(chord.hrdValue()))
                            .brca1Value(chord.brca1Value())
                            .brca2Value(chord.brca2Value())
                            .hrdType(chord.hrdType())
                            .lohCopyNumbers(lohGainDeletions)
                            .drivingGenes(drivingGenes)
                            .build())
                    .build();
        }
        else
        {
            Set<FindingStatus.ResultIssue> errors = new HashSet<>();
            if(experimentType == ExperimentType.TARGETED)
            {
                errors.add(WGS_REQUIRED);
            }
            if(!hasRefSample)
            {
                errors.add(REF_REQUIRED);
            }
            if(!errors.isEmpty())
            {
                return FindingUtil.notAvailableFindingItem(errors);
            }
            else
            {
                throw new IllegalStateException("Chord record should not be null");
            }
        }
    }

    @NotNull
    private static HomologousRecombination.Status hrStatus(@NotNull ChordRecord chord)
    {
        return switch(chord.hrStatus())
        {
            case HR_DEFICIENT -> HomologousRecombination.Status.HR_DEFICIENT;
            case HR_PROFICIENT -> HomologousRecombination.Status.HR_PROFICIENT;
            case CANNOT_BE_DETERMINED -> HomologousRecombination.Status.UNDETERMINED;
            default -> throw new IllegalStateException("Chord status should not be UNKNOWN");
        };
    }

    private static FindingItem<MicrosatelliteStability> createMicrosatelliteStability(PurpleRecord purple,
            LinxRecord linx, DriverFindingList<GainDeletion> gainDeletions, FindingStatus findingStatus)
    {
        MicrosatelliteStability.Status microsatelliteStatus =
                microsatelliteStatus(purple.characteristics().microsatelliteStatus());
        boolean isPresent = microsatelliteStatus == MicrosatelliteStability.Status.MSI;
        List<GainDeletion> lohGainDeletions = isPresent
                ? filterLohGainDeletions(gainDeletions, Genes.MSI_GENES)
                : List.of();
        List<String> drivingGenes = isPresent ? GeneListUtil.genes(purple.reportableSomaticVariants(),
                purple.reportableSomaticGainsDels(),
                linx.germlineHomozygousDisruptions(),
                Genes.MSI_GENES) : List.of();

        return FindingItemBuilder.<MicrosatelliteStability>builder()
                .status(findingStatus)
                .finding(MicrosatelliteStabilityBuilder.builder()
                        .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                        .status(microsatelliteStatus)
                        .indelsPerMb(ThresholdValueFactory.mssValue(purple.characteristics().microsatelliteIndelsPerMb()))
                        .lohCopyNumbers(lohGainDeletions)
                        .drivingGenes(drivingGenes)
                        .build())
                .build();
    }

    @NotNull
    private static MicrosatelliteStability.Status microsatelliteStatus(PurpleMicrosatelliteStatus status)
    {
        return switch(status)
        {
            case MSS -> MicrosatelliteStability.Status.MSS;
            case MSI -> MicrosatelliteStability.Status.MSI;
            case UNKNOWN -> throw new IllegalStateException("Microsatellite status should not be UNKNOWN");
        };
    }

    private static List<GainDeletion> filterLohGainDeletions(
            DriverFindingList<GainDeletion> gainDeletions, Set<String> geneNames)
    {
        return gainDeletions.findings().stream()
                .filter(x -> geneNames.contains(x.gene()))
                .filter(GainDeletion::isLossOfHeterozygosity)
                .sorted(Comparator.comparing(GainDeletion::gene))
                .collect(Collectors.toList());
    }

    public static DriverFindingList<Fusion> createFusionsFindings(LinxRecord linx, FindingStatus findingStatus)
    {
        return DriverFindingListBuilder.<Fusion>builder()
                .status(findingStatus)
                .findings(linx.reportableSomaticFusions().stream()
                        .map(o -> convertFusion(o, DriverSource.SOMATIC)).sorted(Fusion.COMPARATOR).toList())
                .build();
    }

    public static Fusion convertFusion(LinxFusion fusion, DriverSource sampleType)
    {
        DriverInterpretation driverInterpretation = toDriverInterpretation(fusion.driverLikelihood());

        boolean isDriverGene = !fusion.unreportedReasons().contains(LinxUnreportableReason.NOT_KNOWN);
        List<Fusion.UnreportableReason> unreportableReasons = fusion.unreportedReasons().stream()
                .map(o -> Fusion.UnreportableReason.valueOf(o.name()))
                .toList();

        return FusionBuilder.builder()
                .driver(DriverFieldsBuilder.builder()
                        .findingKey(FindingKeys.fusion(sampleType, fusion))
                        .driverSource(sampleType)
                        .reportedStatus(DriverUtil.reportedStatus(isDriverGene, fusion.reported(), driverInterpretation))
                        .driverInterpretation(driverInterpretation)
                        .driverLikelihood(fusion.driverLikelihood() == FusionLikelihoodType.HIGH ? 1.0 : 0.0)
                        .build()
                )
                .geneUp(fusion.geneStart())
                .geneContextUp(fusion.geneContextStart())
                .geneTranscriptUp(fusion.geneTranscriptStart())
                .geneDown(fusion.geneEnd())
                .geneContextDown(fusion.geneContextEnd())
                .geneTranscriptDown(fusion.geneTranscriptEnd())
                .reportedType(Fusion.FusionType.valueOf(fusion.reportedType().name()))
                .unreportedReasons(unreportableReasons)
                .phased(Fusion.FusionPhasedType.valueOf(fusion.phased().name()))
                .fusedExonUp(fusion.fusedExonUp())
                .fusedExonDown(fusion.fusedExonDown())
                .chainLinks(fusion.chainLinks())
                .chainTerminated(fusion.chainTerminated())
                .domainsKept(List.of(fusion.domainsKept().split(";")))
                .domainsLost(List.of(fusion.domainsLost().split(";")))
                .junctionCopyNumber(fusion.junctionCopyNumber())
                .build();
    }

    private static DriverInterpretation toDriverInterpretation(FusionLikelihoodType likelihood)
    {
        return switch(likelihood)
        {
            case HIGH -> DriverInterpretation.HIGH;
            case LOW -> DriverInterpretation.LOW;
            case NA -> DriverInterpretation.UNKNOWN;
        };
    }

    private static DriverFindingList<Virus> createVirusFindings(@Nullable VirusInterpreterData virusInterpreter,
            @NotNull ExperimentType experimentType,
            @NotNull FindingStatus findingStatus)
    {
        if(virusInterpreter != null)
        {
            return DriverFindingListBuilder.<Virus>builder()
                    .status(findingStatus)
                    .findings(convertViruses(virusInterpreter.allViruses()))
                    .build();
        }
        else
        {
            if(experimentType == ExperimentType.TARGETED)
            {
                return FindingUtil.notAvailableDriverFindingList(Set.of(WGS_REQUIRED));
            }
            else
            {
                throw new IllegalStateException("Virus Interpreter data should not be null");
            }
        }
    }

    private static List<Virus> convertViruses(List<VirusInterpreterEntry> viruses)
    {
        return viruses.stream()
                .map(v -> VirusBuilder.builder()
                        .driver(DriverFieldsBuilder.builder()
                                .findingKey(FindingKeys.virus(v))
                                .driverSource(DriverSource.SOMATIC)
                                .reportedStatus(DriverUtil.reportedStatus(
                                        true,
                                        v.reported(),
                                        virusDriverInterpretation(v.driverLikelihood())))
                                .driverInterpretation(virusDriverInterpretation(v.driverLikelihood()))
                                .driverLikelihood(v.driverLikelihood() == VirusLikelihoodType.HIGH ? 1.0 : 0.0)
                                .build()
                        )
                        .name(v.name())
                        .qcStatus(Virus.VirusBreakendQCStatus.valueOf(v.qcStatus().name()))
                        .integrations(v.integrations())
                        .oncogenicVirus(v.interpretation() != null ?
                                Virus.OncogenicVirus.valueOf(Objects.requireNonNull(v.interpretation()).name())
                                : null)
                        .percentageCovered(v.percentageCovered())
                        .meanCoverage(v.meanCoverage())
                        .expectedClonalCoverage(v.expectedClonalCoverage())
                        .build()
                ).sorted(Virus.COMPARATOR).toList();
    }

    private static DriverInterpretation virusDriverInterpretation(VirusLikelihoodType virusLikelihoodType)
    {
        return switch(virusLikelihoodType)
        {
            case LOW -> DriverInterpretation.LOW;
            case HIGH -> DriverInterpretation.HIGH;
            case UNKNOWN -> DriverInterpretation.UNKNOWN;
        };
    }

    private static FindingList<PharmacoGenotype> createPharmacoGenotypesFindings(@Nullable Set<PeachGenotype> peachGenotypes,
            @NotNull FindingStatus findingStatus)
    {
        if(peachGenotypes != null)
        {
            return FindingListBuilder.<PharmacoGenotype>builder()
                    .status(findingStatus)
                    .findings(peachGenotypes.stream().map(o ->
                                    PharmacoGenotypeBuilder.builder()
                                            .findingKey(FindingKeys.pharmacoGenotype(o.gene(), o.allele()))
                                            .gene(o.gene())
                                            .allele(o.allele())
                                            .alleleCount(o.alleleCount())
                                            .function(o.function())
                                            .haplotype(o.haplotype())
                                            .linkedDrugs(o.linkedDrugs())
                                            .urlPrescriptionInfo(o.urlPrescriptionInfo())
                                            .build())
                            .sorted(PharmacoGenotype.COMPARATOR)
                            .toList())
                    .build();
        }
        else
        {
            throw new IllegalStateException("Peach genotypes should not be null");
        }
    }
}
