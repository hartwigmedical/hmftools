package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.DisruptionFactory.createGermlineDisruptions;
import static com.hartwig.hmftools.finding.DisruptionFactory.createSomaticDisruptions;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
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
import com.hartwig.hmftools.finding.datamodel.DriverFieldsBuilder;
import com.hartwig.hmftools.finding.datamodel.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.DriverInterpretation;
import com.hartwig.hmftools.finding.datamodel.DriverSource;
import com.hartwig.hmftools.finding.datamodel.FindingItem;
import com.hartwig.hmftools.finding.datamodel.FindingItemBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingList;
import com.hartwig.hmftools.finding.datamodel.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.FindingsStatus;
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
import com.hartwig.hmftools.finding.datamodel.VisualisationFiles;

import org.jetbrains.annotations.Nullable;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory
{
    public static FindingRecord fromOrangeJsonWithTranscriptFile(Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        try(Reader reader = Files.newBufferedReader(orangeJson))
        {
            OrangeRecord orangeRecord = OrangeJson.getInstance().read(reader);
            return fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv);
        }
    }

    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        FindingConfig findingConfig =
                FindingConfig.createFindingConfig(clinicalTranscriptsTsv, driverGeneTsv, orangeRecord.refGenomeVersion());

        LinxRecord linx = orangeRecord.linx();
        PurpleRecord purple = orangeRecord.purple();

        Set<PurpleQCStatus> purpleQCStatuses = purple.fit().qc().status();
        boolean hasContamination = purpleQCStatuses.contains(PurpleQCStatus.FAIL_CONTAMINATION);
        boolean hasReliablePurity =
                !purpleQCStatuses.contains(PurpleQCStatus.FAIL_NO_TUMOR) && !purpleQCStatuses.contains(PurpleQCStatus.WARN_LOW_PURITY);

        FindingRecordBuilder
                builder = FindingRecordBuilder.builder()
                .version("1.0")
                .metaProperties(MetaPropertiesBuilder.builder()
                        .refGenomeVersion(RefGenomeVersion.valueOf(orangeRecord.refGenomeVersion().name()))
                        .sequencingScope(SequencingScope.valueOf(orangeRecord.experimentType().name()))
                        .pipelineVersion(orangeRecord.pipelineVersion())
                        .sampleId(orangeRecord.sampleId())
                        .samplingDate(orangeRecord.samplingDate())
                        .build())
                .qc(createQc(purple))
                .fusions(createFusionsFindings(orangeRecord.linx()));

        DriverFindingList<GainDeletion>
                somaticGainDeletions = addPurpleFindings(builder, orangeRecord, findingConfig, hasReliablePurity);

        if(hasReliablePurity)
        {
            builder.somaticDisruptions(createSomaticDisruptions(linx))
                    .germlineDisruptions(createGermlineDisruptions(orangeRecord.refSample() != null, linx))
                    .viruses(createVirusFindings(orangeRecord.virusInterpreter()))
                    .homologousRecombination(createHomologousRecombination(orangeRecord.chord(), purple, linx, somaticGainDeletions));
        }
        else
        {
            builder.somaticDisruptions(FindingUtil.notAvailableDriverFindingList())
                    .germlineDisruptions(FindingUtil.notAvailableDriverFindingList())
                    .viruses(FindingUtil.notAvailableDriverFindingList())
                    .homologousRecombination(FindingUtil.notAvailableFindingItem());
        }

        VisualisationFiles visualisationFiles = VisualisationFilesFactory.create(orangeRecord.plots());

        return builder.predictedTumorOrigins(createPredictedTumorOriginList(orangeRecord.cuppa()))
                .hlaAlleles(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, hasReliablePurity, hasContamination))
                .pharmacoGenotypes(createPharmacoGenotypesFindings(orangeRecord.peach(), hasContamination))
                .visualisationFiles(visualisationFiles)
                .build();
    }

    // return the gain deletions cause they are needed by HRD, will see if we can find a better way
    private static DriverFindingList<GainDeletion> addPurpleFindings(
            FindingRecordBuilder builder, final OrangeRecord orangeRecord,
            FindingConfig findingConfig,
            boolean hasReliablePurity)
    {
        boolean hasRefSample = orangeRecord.refSample() != null;

        PurpleRecord purple = orangeRecord.purple();

        builder.purityPloidyFit(createPurityPloidyFit(purple));

        DriverFindingList<GainDeletion> somaticGainDeletions;

        DriverFindingList<SmallVariant> smallVariants =
                SmallVariantFactory.somaticSmallVariantFindings(purple, FindingsStatus.OK, findingConfig);
        if(hasReliablePurity)
        {
            somaticGainDeletions =
                    GainDeletionFactory.somaticGainDeletionFindings(orangeRecord.refGenomeVersion(), FindingsStatus.OK, purple);

            builder.somaticSmallVariants(smallVariants)
                    .germlineSmallVariants(SmallVariantFactory.germlineSmallVariantFindings(hasRefSample, purple, findingConfig))
                    .somaticGainDeletions(somaticGainDeletions)
                    .germlineGainDeletions(GainDeletionFactory.germlineGainDeletionFindings(hasRefSample, orangeRecord.refGenomeVersion(), purple))
                    .microsatelliteStability(createMicrosatelliteStability(purple, orangeRecord.linx(), somaticGainDeletions))
                    .tumorMutationalLoad(createTumorMutationalLoad(purple))
                    .tumorMutationalBurden(createTumorMutationalBurden(purple));

        }
        else
        {
            somaticGainDeletions = FindingUtil.notAvailableDriverFindingList();
            builder.somaticSmallVariants(smallVariants)
                    .germlineSmallVariants(FindingUtil.notAvailableDriverFindingList())
                    .somaticGainDeletions(somaticGainDeletions)
                    .germlineGainDeletions(FindingUtil.notAvailableDriverFindingList())
                    .microsatelliteStability(FindingUtil.notAvailableFindingItem())
                    .tumorMutationalLoad(FindingUtil.notAvailableFindingItem())
                    .tumorMutationalBurden(FindingUtil.notAvailableFindingItem());
        }
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

    private static PurityPloidyFit createPurityPloidyFit(PurpleRecord purple)
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
                .build();
    }

    private static FindingList<PredictedTumorOrigin> createPredictedTumorOriginList(@Nullable CuppaData cuppa)
    {
        if(cuppa != null)
        {
            return FindingListBuilder.<PredictedTumorOrigin>builder()
                    .status(FindingsStatus.OK)
                    .findings(List.of(createPredictedTumorOrigin(cuppa.bestPrediction(), cuppaMode(cuppa.mode()))))
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingList();
        }
    }

    private static PredictedTumorOrigin createPredictedTumorOrigin(CuppaPrediction prediction, PredictedTumorOrigin.CuppaMode cuppaMode)
    {
        return PredictedTumorOriginBuilder.builder()
                .findingKey(FindingKeys.predictedTumorOrigin(prediction.cancerType()))
                .mode(cuppaMode)
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

    private static FindingItem<TumorMutationalBurden> createTumorMutationalBurden(PurpleRecord purple)
    {
        TumorMutationalBurden.Status status = tumorMutationalBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus());
        if(status != null)
        {
            return FindingItemBuilder.<TumorMutationalBurden>builder()
                    .status(FindingsStatus.OK)
                    .finding(TumorMutationalBurdenBuilder.builder()
                            .findingKey(FindingKeys.tumorMutationBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus()))
                            .burdenPerMb(purple.characteristics().tumorMutationalBurdenPerMb())
                            .status(status)
                            .svBurden(purple.characteristics().svTumorMutationalBurden())
                            .build())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingItem();
        }
    }

    @Nullable
    private static TumorMutationalBurden.Status tumorMutationalBurdenStatus(PurpleTumorMutationalStatus status)
    {
        return switch(status)
        {
            case HIGH -> TumorMutationalBurden.Status.HIGH;
            case LOW -> TumorMutationalBurden.Status.LOW;
            case UNKNOWN -> null;
        };
    }

    private static FindingItem<TumorMutationalLoad> createTumorMutationalLoad(PurpleRecord purple)
    {
        TumorMutationalLoad.Status status = tumorMutationalLoadStatus(purple.characteristics().tumorMutationalLoadStatus());
        if(status != null)
        {
            return FindingItemBuilder.<TumorMutationalLoad>builder()
                    .status(FindingsStatus.OK)
                    .finding(TumorMutationalLoadBuilder.builder()
                            .findingKey(FindingKeys.tumorMutationLoadStatus(purple.characteristics().tumorMutationalLoadStatus()))
                            .load(purple.characteristics().tumorMutationalLoad())
                            .status(status)
                            .build())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingItem();
        }
    }

    @Nullable
    private static TumorMutationalLoad.Status tumorMutationalLoadStatus(PurpleTumorMutationalStatus status)
    {
        return switch(status)
        {
            case HIGH -> TumorMutationalLoad.Status.HIGH;
            case LOW -> TumorMutationalLoad.Status.LOW;
            case UNKNOWN -> null;
        };
    }

    private static FindingItem<HomologousRecombination> createHomologousRecombination(@Nullable ChordRecord chord,
            PurpleRecord purple,
            LinxRecord linx,
            DriverFindingList<GainDeletion> gainDeletions)
    {
        // TODO: Should distinguish between unknown and cannot be determined?
        HomologousRecombination.HrStatus hrStatus = hrStatus(chord);
        if(hrStatus != null)
        {
            List<GainDeletion> lohGainDeletions = hrStatus == HomologousRecombination.HrStatus.HR_DEFICIENT
                    ? filterLohGainDeletions(gainDeletions, Genes.HRD_GENES)
                    : List.of();
            return FindingItemBuilder.<HomologousRecombination>builder()
                    .status(FindingsStatus.OK)
                    .finding(HomologousRecombinationBuilder.builder()
                            .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                            .brca1Value(chord.brca1Value())
                            .brca2Value(chord.brca2Value())
                            .hrdValue(chord.hrdValue())
                            .hrStatus(hrStatus)
                            .hrdType(chord.hrdType())
                            .lohCopyNumbers(lohGainDeletions)
                            .relatedGenes(GeneListUtil.genes(purple.reportableSomaticVariants(),
                                    purple.reportableSomaticGainsDels(),
                                    linx.germlineHomozygousDisruptions(),
                                    Genes.HRD_GENES).stream().toList())
                            .build())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingItem();
        }
    }

    @Nullable
    private static HomologousRecombination.HrStatus hrStatus(@Nullable ChordRecord chord)
    {
        return chord != null ? switch(chord.hrStatus())
        {
            case HR_DEFICIENT -> HomologousRecombination.HrStatus.HR_DEFICIENT;
            case HR_PROFICIENT -> HomologousRecombination.HrStatus.HR_PROFICIENT;
            default -> null;
        } : null;
    }

    private static FindingItem<MicrosatelliteStability> createMicrosatelliteStability(PurpleRecord purple,
            LinxRecord linx, DriverFindingList<GainDeletion> gainDeletions)
    {
        MicrosatelliteStability.MicrosatelliteStatus microsatelliteStatus =
                microsatelliteStatus(purple.characteristics().microsatelliteStatus());
        if(microsatelliteStatus != null)
        {
            List<GainDeletion> lohGainDeletions = microsatelliteStatus == MicrosatelliteStability.MicrosatelliteStatus.MSI
                    ? filterLohGainDeletions(gainDeletions, Genes.MSI_GENES)
                    : List.of();
            return FindingItemBuilder.<MicrosatelliteStability>builder()
                    .status(FindingsStatus.OK)
                    .finding(MicrosatelliteStabilityBuilder.builder()
                            .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                            .microsatelliteStatus(microsatelliteStatus)
                            .microsatelliteIndelsPerMb(purple.characteristics().microsatelliteIndelsPerMb())
                            .lohCopyNumbers(lohGainDeletions)
                            .relatedGenes(GeneListUtil.genes(purple.reportableSomaticVariants(),
                                    purple.reportableSomaticGainsDels(),
                                    linx.germlineHomozygousDisruptions(),
                                    Genes.MSI_GENES).stream().toList())
                            .build())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingItem();
        }
    }

    @Nullable
    private static MicrosatelliteStability.MicrosatelliteStatus microsatelliteStatus(PurpleMicrosatelliteStatus status)
    {
        return switch(status)
        {
            case MSS -> MicrosatelliteStability.MicrosatelliteStatus.MSS;
            case MSI -> MicrosatelliteStability.MicrosatelliteStatus.MSI;
            default -> null;
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

    public static DriverFindingList<Fusion> createFusionsFindings(LinxRecord linx)
    {
        return DriverFindingListBuilder.<Fusion>builder()
                .status(FindingsStatus.OK)
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
                .geneStart(fusion.geneStart())
                .geneContextStart(fusion.geneContextStart())
                .geneTranscriptStart(fusion.geneTranscriptStart())
                .geneEnd(fusion.geneEnd())
                .geneContextEnd(fusion.geneContextEnd())
                .geneTranscriptEnd(fusion.geneTranscriptEnd())
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

    private static DriverFindingList<Virus> createVirusFindings(@Nullable VirusInterpreterData virusInterpreter)
    {
        if(virusInterpreter != null)
        {
            return DriverFindingListBuilder.<Virus>builder()
                    .status(FindingsStatus.OK)
                    .findings(convertViruses(virusInterpreter.allViruses()))
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableDriverFindingList();
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
            boolean hasContamination)
    {
        if(peachGenotypes != null)
        {
            return FindingListBuilder.<PharmacoGenotype>builder()
                    .status(hasContamination ? FindingsStatus.NOT_RELIABLE : FindingsStatus.OK)
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
            return FindingUtil.notAvailableFindingList();
        }
    }
}
