package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.DisruptionFactory.createGermlineDisruptions;
import static com.hartwig.hmftools.finding.DisruptionFactory.createSomaticDisruptions;
import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.NORMAL_REQUIRED;
import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.WGS_REQUIRED;
import static com.hartwig.hmftools.finding.util.FindingUtil.germlineStatus;
import static com.hartwig.hmftools.finding.util.FindingUtil.notAvailableDriverFindingList;
import static com.hartwig.hmftools.finding.util.FindingUtil.notAvailableFindingItem;
import static com.hartwig.hmftools.finding.util.FindingUtil.notAvailableFindingList;
import static com.hartwig.hmftools.finding.util.FindingUtil.somaticStatus;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
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
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.Fusion;
import com.hartwig.hmftools.finding.datamodel.FusionBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombination;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombinationBuilder;
import com.hartwig.hmftools.finding.datamodel.HomologousRecombinationPredictionBuilder;
import com.hartwig.hmftools.finding.datamodel.MetaProperties;
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
import com.hartwig.hmftools.finding.datamodel.driver.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.driver.DriverSource;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import org.jetbrains.annotations.Nullable;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory
{
    public static FindingRecord fromOrangeJsonWithTranscriptFile(Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, @Nullable Gender gender) throws IOException
    {
        try(Reader reader = Files.newBufferedReader(orangeJson))
        {
            OrangeRecord orangeRecord = OrangeJson.getInstance().read(reader);
            return fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv, gender);
        }
    }

    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv, @Nullable Gender gender) throws IOException
    {
        FindingConfig findingConfig =
                FindingConfig.createFindingConfig(clinicalTranscriptsTsv, driverGeneTsv, orangeRecord.refGenomeVersion(), gender, true);

        LinxRecord linx = orangeRecord.linx();
        PurpleRecord purple = orangeRecord.purple();

        Qc qc = createQc(purple);

        FindingStatus findingStatus = FindingsStatusFactory.toFindingsStatus(qc.errors(), qc.warnings());

        boolean hasRefSample = orangeRecord.refSample() != null;

        ExperimentType experimentType = orangeRecord.experimentType();

        DriverFindingList<SmallVariant> smallVariants =
                SmallVariantFactory.somaticSmallVariantFindings(purple, findingStatus, findingConfig);

        ArmCopyNumberFactory cnPerChromosome = new ArmCopyNumberFactory(
                purple.allSomaticCopyNumbers(), purple.fit().ploidy(), findingConfig.gender(), orangeRecord.refGenomeVersion());

        DriverFindingList<GainDeletion> somaticGainDeletions =
                GainDeletionFactory.somaticGainDeletionFindings(findingStatus, purple, cnPerChromosome, findingConfig.geneCopyNumbersOptional());

        DriverFindingList<Disruption> germlineDisruptions =
                createGermlineDisruptions(orangeRecord.refSample() != null, linx, findingStatus);

        List<Disruption> germlineHomozygousDisruptions =
                germlineDisruptions.findingsIfOk().stream().filter(Disruption::isHomozygous).toList();

        MetaProperties metaProperties = createMetaProperties(orangeRecord, experimentType);

        String version = FindingRecord.VERSION != null ? FindingRecord.VERSION : "local-dev";
        return FindingRecordBuilder.builder()
                .version(version)
                .metaProperties(metaProperties)
                .qc(qc)
                .purityPloidyFit(createPurityPloidyFit(purple, orangeRecord.plots()))
                .fusions(createFusionsFindings(orangeRecord.linx(), findingStatus))
                .somaticSmallVariants(smallVariants)
                .germlineSmallVariants(SmallVariantFactory.germlineSmallVariantFindings(hasRefSample, purple, findingStatus, findingConfig))
                .somaticGainDeletions(somaticGainDeletions)
                .germlineGainDeletions(GainDeletionFactory.germlineGainDeletionFindings(hasRefSample, findingStatus, purple, cnPerChromosome, findingConfig.geneCopyNumbersOptional()))
                .microsatelliteStability(createMicrosatelliteStability(purple,
                        smallVariants,
                        somaticGainDeletions,
                        germlineHomozygousDisruptions,
                        findingStatus
                ))
                .tumorMutationalLoad(createTumorMutationalLoad(purple, findingStatus))
                .tumorMutationalBurden(createTumorMutationalBurden(purple, findingStatus))
                .chromosomeArmCopyNumbers(cnPerChromosome.toArmCopyNumberFindings(findingStatus))
                .somaticDisruptions(createSomaticDisruptions(linx, findingStatus))
                .germlineDisruptions(germlineDisruptions)
                .viruses(createVirusFindings(orangeRecord.virusInterpreter(), experimentType, findingStatus))
                .homologousRecombination(createHomologousRecombination(orangeRecord.chord(),
                        smallVariants,
                        somaticGainDeletions,
                        germlineHomozygousDisruptions,
                        findingStatus, experimentType,
                        hasRefSample
                ))
                .predictedTumorOrigin(createPredictedTumorOrigin(orangeRecord.cuppa(), orangeRecord.plots(), experimentType, findingStatus))
                .hlaAlleles(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, findingStatus))
                .pharmacoGenotypes(createPharmacoGenotypesFindings(orangeRecord.peach(), findingStatus))
                .highExpressionGenes(notAvailableFindingList(Set.of(FindingStatus.Issue.RNA_REQUIRED)))
                .lowExpressionGenes(notAvailableFindingList(Set.of(FindingStatus.Issue.RNA_REQUIRED)))
                .rnaFusions(notAvailableFindingList(Set.of(FindingStatus.Issue.RNA_REQUIRED)))
                .novelSpliceJunctions(notAvailableFindingList(Set.of(FindingStatus.Issue.RNA_REQUIRED)))
                .build();
    }

    private static MetaProperties createMetaProperties(final OrangeRecord orangeRecord, final ExperimentType experimentType)
    {
        return MetaPropertiesBuilder.builder()
                .refGenomeVersion(RefGenomeVersion.valueOf(orangeRecord.refGenomeVersion().name()))
                .sequencingScope(SequencingScope.valueOf(experimentType.name()))
                .pipelineVersion(orangeRecord.pipelineVersion())
                .sampleId(orangeRecord.sampleId())
                .samplingDate(orangeRecord.samplingDate())
                .build();
    }

    private static Qc createQc(PurpleRecord purple)
    {
        PurpleFit purpleFit = purple.fit();
        SortedSet<Qc.GermlineAberration> germlineAberrations = purpleFit.qc().germlineAberrations().stream()
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
                .collect(Collectors.toCollection(TreeSet::new));

        return QcBuilder.builder()
                .isPass(purpleFit.qc().status().contains(PurpleQCStatus.PASS))
                .errors(QcStatusFactory.toErrors(purpleFit.qc().status()))
                .warnings(QcStatusFactory.toWarnings(purpleFit.qc().status()))
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
                .purpleRainfallPlot(VisualisationFileUtil.create(orangePlots.purpleKataegisPlot()))
                .build();
    }

    private static FindingItem<PredictedTumorOrigin> createPredictedTumorOrigin(@Nullable CuppaData cuppa, OrangePlots orangePlots,
            ExperimentType experimentType,
            FindingStatus findingStatus)
    {
        if(cuppa != null)
        {
            List<PredictedTumorOrigin.Prediction> predictedTumorOrigins = cuppa.predictions().stream()
                    .map(FindingRecordFactory::createPredictedTumorOriginPrediction)
                    .toList();
            return FindingItemBuilder.<PredictedTumorOrigin>builder()
                    .status(findingStatus)
                    .finding(PredictedTumorOriginBuilder.builder()
                            .findingKey("predictedTumorOrigin")
                            .mode(cuppaMode(cuppa.mode()))
                            .predictions(predictedTumorOrigins)
                            .visualisationFile(VisualisationFileUtil.createNullable(orangePlots.cuppaSummaryPlot()))
                            .build()
                    )
                    .build();
        }
        else
        {
            if(experimentType == ExperimentType.TARGETED)
            {
                return notAvailableFindingItem(Set.of(WGS_REQUIRED));
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
                .status(somaticStatus(findingStatus))
                .finding(TumorMutationalLoadBuilder.builder()
                        .findingKey(FindingKeys.tumorMutationLoadStatus(purple.characteristics().tumorMutationalLoadStatus()))
                        .load(ThresholdValueFactory.tmlValue(purple.characteristics().tumorMutationalLoad()))
                        .status(status)
                        .build())
                .build();
    }

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
            DriverFindingList<SmallVariant> smallVariants,
            DriverFindingList<GainDeletion> gainDeletions,
            List<Disruption> germlineHomozygousDisruptions,
            FindingStatus findingStatus,
            ExperimentType experimentType,
            boolean hasRefSample)
    {
        if(chord != null)
        {
            HomologousRecombination.Prediction prediction = createHomologousRecombinationPrediction(chord);
            boolean isPresent = prediction.status() == HomologousRecombination.Status.HR_DEFICIENT;
            SortedSet<String> drivingGenes = isPresent ? GeneListUtil.genes(smallVariants,
                    gainDeletions,
                    germlineHomozygousDisruptions,
                    Genes.HRD_GENES) : new TreeSet<>();
            SortedMap<HomologousRecombination.HrdCancerType, HomologousRecombination.Prediction> predictions = new TreeMap<>();
            predictions.put(prediction.cancerType(), prediction);
            return FindingItemBuilder.<HomologousRecombination>builder()
                    .status(findingStatus)
                    .finding(HomologousRecombinationBuilder.builder()
                            .drivingGenes(drivingGenes)
                            .predictions(predictions)
                            .build())
                    .build();
        }
        else
        {
            Set<FindingStatus.Issue> errors = new HashSet<>();
            if(experimentType == ExperimentType.TARGETED)
            {
                errors.add(WGS_REQUIRED);
            }
            if(!hasRefSample)
            {
                errors.add(NORMAL_REQUIRED);
            }
            if(!errors.isEmpty())
            {
                return notAvailableFindingItem(errors);
            }
            else
            {
                throw new IllegalStateException("Chord record should not be null");
            }
        }
    }

    private static HomologousRecombination.Prediction createHomologousRecombinationPrediction(ChordRecord chord)
    {
        HomologousRecombination.Status hrStatus = hrStatus(chord);
        HomologousRecombination.HrdType hrdType = switch (chord.hrdType())
        {
            case "none" -> null;
            case "BRCA1_type" -> HomologousRecombination.HrdType.BRACA1_TYPE;
            case "BRCA2_type" -> HomologousRecombination.HrdType.BRACA2_TYPE;
            case "cannot_be_determined" -> HomologousRecombination.HrdType.CANNOT_BE_DETERMINED;
            default -> throw new IllegalStateException("Unexpected chord hrdType: " + chord.hrdType());
        };
        return HomologousRecombinationPredictionBuilder.builder()
                .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                .status(hrStatus)
                .hrdProbability(ThresholdValueFactory.hrdValue(chord.hrdValue()))
                .brca1Probability(chord.brca1Value())
                .brca2Probability(chord.brca2Value())
                .hrdType(hrdType)
                .cancerType(HomologousRecombination.HrdCancerType.PAN_CANCER)
                .build();
    }



    private static HomologousRecombination.Status hrStatus(ChordRecord chord)
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
            DriverFindingList<SmallVariant> smallVariants, DriverFindingList<GainDeletion> gainDeletions,
            List<Disruption> germlineHomozygousDisruptions,
            FindingStatus findingStatus)
    {
        MicrosatelliteStability.Status microsatelliteStatus =
                microsatelliteStatus(purple.characteristics().microsatelliteStatus());
        boolean isPresent = microsatelliteStatus == MicrosatelliteStability.Status.MSI;
        SortedSet<String> drivingGenes = isPresent ? GeneListUtil.genes(smallVariants,
                gainDeletions,
                germlineHomozygousDisruptions,
                Genes.MSI_GENES) : new TreeSet<>();

        return FindingItemBuilder.<MicrosatelliteStability>builder()
                .status(somaticStatus(findingStatus))
                .finding(MicrosatelliteStabilityBuilder.builder()
                        .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                        .status(microsatelliteStatus)
                        .indelsPerMb(ThresholdValueFactory.msiValue(purple.characteristics().microsatelliteIndelsPerMb()))
                        .drivingGenes(drivingGenes)
                        .build())
                .build();
    }

    private static MicrosatelliteStability.Status microsatelliteStatus(PurpleMicrosatelliteStatus status)
    {
        return switch(status)
        {
            case MSS -> MicrosatelliteStability.Status.MSS;
            case MSI -> MicrosatelliteStability.Status.MSI;
            case UNKNOWN -> throw new IllegalStateException("Microsatellite status should not be UNKNOWN");
        };
    }

    public static DriverFindingList<Fusion> createFusionsFindings(LinxRecord linx, FindingStatus findingStatus)
    {
        return DriverFindingListBuilder.<Fusion>builder()
                .status(somaticStatus(findingStatus))
                .findings(linx.reportableSomaticFusions().stream()
                        .map(o -> convertFusion(o, DriverSource.SOMATIC)).sorted(Fusion.COMPARATOR).toList())
                .build();
    }

    public static Fusion convertFusion(LinxFusion fusion, DriverSource sampleType)
    {
        DriverInterpretation driverInterpretation = toDriverInterpretation(fusion.driverLikelihood());

        boolean isDriverGene = !fusion.unreportedReasons().contains(LinxUnreportableReason.NOT_KNOWN);

        Fusion.FusionType reportedType = switch(fusion.reportedType())
        {
            case NONE -> Fusion.FusionType.NONE;
            case PROMISCUOUS_3 -> Fusion.FusionType.PROMISCUOUS_3;
            case PROMISCUOUS_5 -> Fusion.FusionType.PROMISCUOUS_5;
            case PROMISCUOUS_BOTH -> Fusion.FusionType.PROMISCUOUS_BOTH;
            case IG_PROMISCUOUS -> Fusion.FusionType.ENHANCER_PROMISCUOUS;
            case KNOWN_PAIR -> Fusion.FusionType.KNOWN_PAIR;
            case IG_KNOWN_PAIR -> Fusion.FusionType.ENHANCER_KNOWN_PAIR;
            case EXON_DEL_DUP -> Fusion.FusionType.EXON_DEL_DUP;
            case PROMISCUOUS_ENHANCER_TARGET -> Fusion.FusionType.PROMISCUOUS_ENHANCER_TARGET;
        };

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
                .reportedType(reportedType)
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
            ExperimentType experimentType,
            FindingStatus findingStatus)
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
                return notAvailableDriverFindingList(Set.of(WGS_REQUIRED));
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
                        .percentageCovered(toPercentage(v.percentageCovered()))
                        .meanCoverage(v.meanCoverage())
                        .expectedClonalCoverage(toPercentageNullable(v.expectedClonalCoverage()))
                        .build()
                ).sorted(Virus.COMPARATOR).toList();
    }

    @Nullable
    private static Double toPercentageNullable(@Nullable Double value)
    {
        return value == null ? null : toPercentage(value);
    }

    private static double toPercentage(double value)
    {
        return value * 0.01;
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
            FindingStatus findingStatus)
    {
        if(peachGenotypes != null)
        {
            return FindingListBuilder.<PharmacoGenotype>builder()
                    .status(germlineStatus(findingStatus))
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
