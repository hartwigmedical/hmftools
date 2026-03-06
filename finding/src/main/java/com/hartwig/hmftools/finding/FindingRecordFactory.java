package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.DisruptionFactory.createGermlineDisruptions;
import static com.hartwig.hmftools.finding.DisruptionFactory.createSomaticDisruptions;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
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
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;
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
import com.hartwig.hmftools.finding.datamodel.PharmocoGenotype;
import com.hartwig.hmftools.finding.datamodel.PharmocoGenotypeBuilder;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOrigin;
import com.hartwig.hmftools.finding.datamodel.PredictedTumorOriginBuilder;
import com.hartwig.hmftools.finding.datamodel.PurityPloidyFit;
import com.hartwig.hmftools.finding.datamodel.PurityPloidyFitBuilder;
import com.hartwig.hmftools.finding.datamodel.PurityPloidyFitQcBuilder;
import com.hartwig.hmftools.finding.datamodel.RefGenomeVersion;
import com.hartwig.hmftools.finding.datamodel.SequencingScope;
import com.hartwig.hmftools.finding.datamodel.TumorMutationStatus;
import com.hartwig.hmftools.finding.datamodel.TumorMutationStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.Virus;
import com.hartwig.hmftools.finding.datamodel.VirusBuilder;
import com.hartwig.hmftools.finding.datamodel.VisualisationFiles;

import org.jetbrains.annotations.NotNull;
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
            OrangeRecord orangeRecord = com.hartwig.hmftools.datamodel.OrangeJson.getInstance().read(reader);
            return fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv);
        }
    }

    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRecord.refGenomeVersion(), clinicalTranscriptsTsv) : null;
        Map<String, DriverGene> driverGenes = driverGenesMap(driverGeneTsv);

        LinxRecord linx = orangeRecord.linx();
        PurpleRecord purple = orangeRecord.purple();

        boolean hasReliablePurity = purple.fit().containsTumorCells();
        boolean hasContamination = purple.fit().qc().status().contains(PurpleQCStatus.FAIL_CONTAMINATION);

        FindingRecordBuilder
                builder = FindingRecordBuilder.builder()
                .metaProperties(MetaPropertiesBuilder.builder()
                        .version("1.0")
                        .refGenomeVersion(RefGenomeVersion.valueOf(orangeRecord.refGenomeVersion().name()))
                        .sequencingScope(SequencingScope.valueOf(orangeRecord.experimentType().name()))
                        .pipelineVersion(orangeRecord.pipelineVersion())
                        .sampleId(orangeRecord.sampleId())
                        .samplingDate(orangeRecord.samplingDate())
                        .build())
                .somaticDisruptions(createSomaticDisruptions(hasReliablePurity, linx))
                .germlineDisruptions(createGermlineDisruptions(orangeRecord.refSample() != null, linx))
                .fusions(createFusionsFindings(orangeRecord.linx()));

        DriverFindingList<GainDeletion>
                somaticGainDeletions = addPurpleFindings(builder, orangeRecord, clinicalTranscriptsModel, driverGenes);

        VisualisationFiles visualisationFiles = VisualisationFilesFactory.create(orangeRecord.plots());

        return builder.predictedTumorOrigins(createPredictedTumorOriginList(orangeRecord.cuppa()))
                .homologousRecombination(createHomologousRecombination(orangeRecord.chord(), purple, linx, somaticGainDeletions))
                .viruses(createVirusFindings(orangeRecord.virusInterpreter()))
                .hlaAlleles(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, hasReliablePurity))
                .pharmocoGenotypes(createPharmcoGenotypesFindings(orangeRecord.peach(), hasContamination))
                .visualisationFiles(visualisationFiles)
                .build();
    }

    // return the gain deletions cause they are needed by HRD, will see if we can find a better way
    private static DriverFindingList<GainDeletion> addPurpleFindings(
            FindingRecordBuilder builder, final OrangeRecord orangeRecord,
            final @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel, Map<String, DriverGene> driverGenes)
    {
        boolean hasRefSample = orangeRecord.refSample() != null;

        PurpleRecord purple = orangeRecord.purple();

        builder.purityPloidyFit(createPurityPloidyFit(purple));

        DriverFindingList<GainDeletion> somaticGainDeletions;
        FindingsStatus findingsStatus = purpleFindingsStatus(purple);
        if(findingsStatus == FindingsStatus.OK)
        {
            somaticGainDeletions =
                    GainDeletionFactory.somaticGainDeletionFindings(orangeRecord.refGenomeVersion(), findingsStatus, purple);

            builder.somaticSmallVariants(SmallVariantFactory.somaticSmallVariantFindings(purple, findingsStatus, clinicalTranscriptsModel, driverGenes))
                    .germlineSmallVariants(SmallVariantFactory.germlineSmallVariantFindings(hasRefSample, purple, clinicalTranscriptsModel, driverGenes))
                    .somaticGainDeletions(somaticGainDeletions)
                    .germlineGainDeletions(GainDeletionFactory.germlineGainDeletionFindings(hasRefSample, orangeRecord.refGenomeVersion(), purple))
                    .microsatelliteStability(createMicrosatelliteStability(purple, orangeRecord.linx(), somaticGainDeletions))
                    .tumorMutationStatus(createTumorMutationStatus(purple));
        }
        else
        {
            somaticGainDeletions = FindingUtil.notAvailableDriverFindingList();
            builder.somaticSmallVariants(FindingUtil.notAvailableDriverFindingList())
                    .germlineSmallVariants(FindingUtil.notAvailableDriverFindingList())
                    .somaticGainDeletions(somaticGainDeletions)
                    .germlineGainDeletions(FindingUtil.notAvailableDriverFindingList())
                    .microsatelliteStability(FindingUtil.notAvailableFindingItem())
                    .tumorMutationStatus(FindingUtil.notAvailableFindingItem());
        }
        return somaticGainDeletions;
    }

    @NotNull
    private static PurityPloidyFit createPurityPloidyFit(@NotNull PurpleRecord purple)
    {
        PurpleFit purpleFit = purple.fit();
        Set<PurityPloidyFit.QCStatus> qcStatuses = purpleFit.qc().status().stream()
                .map(o -> PurityPloidyFit.QCStatus.valueOf(o.name()))
                .collect(Collectors.toSet());
        Set<PurityPloidyFit.GermlineAberration> germlineAberrations = purpleFit.qc().germlineAberrations().stream()
                .map(o -> PurityPloidyFit.GermlineAberration.valueOf(o.name()))
                .collect(Collectors.toSet());

        return PurityPloidyFitBuilder.builder()
                .qc(PurityPloidyFitQcBuilder.builder()
                        .status(qcStatuses)
                        .germlineAberrations(germlineAberrations)
                        .amberMeanDepth(purpleFit.qc().amberMeanDepth())
                        .contamination(purpleFit.qc().contamination())
                        .totalCopyNumberSegments(purpleFit.qc().totalCopyNumberSegments())
                        .unsupportedCopyNumberSegments(purpleFit.qc().unsupportedCopyNumberSegments())
                        .deletedGenes(purpleFit.qc().deletedGenes())
                        .build())
                .fittedPurityMethod(PurityPloidyFit.FittedPurityMethod.valueOf(purpleFit.fittedPurityMethod().name()))
                .purity(purpleFit.purity())
                .minPurity(purpleFit.minPurity())
                .maxPurity(purpleFit.maxPurity())
                .ploidy(purpleFit.ploidy())
                .minPloidy(purpleFit.minPloidy())
                .maxPloidy(purpleFit.maxPloidy())
                .build();
    }

    private static FindingItem<TumorMutationStatus> createTumorMutationStatus(PurpleRecord purple)
    {
        return FindingItemBuilder.<TumorMutationStatus>builder()
                .status(FindingsStatus.OK)
                .finding(TumorMutationStatusBuilder.builder()
                        .findingKey(FindingKeys.tumorMutationStatus(purple.characteristics().tumorMutationalBurdenStatus(),
                                purple.characteristics().tumorMutationalLoadStatus()))
                        .tumorMutationalBurdenPerMb(purple.characteristics().tumorMutationalBurdenPerMb())
                        .tumorMutationalBurdenStatus(
                                TumorMutationStatus.Status.valueOf(purple.characteristics().tumorMutationalBurdenStatus().name()))
                        .tumorMutationalLoad(purple.characteristics().tumorMutationalLoad())
                        .tumorMutationalLoadStatus(
                                TumorMutationStatus.Status.valueOf(purple.characteristics().tumorMutationalLoadStatus().name()))
                        .svTumorMutationalBurden(purple.characteristics().svTumorMutationalBurden())
                        .build())
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

    private static FindingItem<HomologousRecombination> createHomologousRecombination(@Nullable ChordRecord chord,
            PurpleRecord purple,
            LinxRecord linx,
            DriverFindingList<GainDeletion> gainDeletions)
    {
        if(chord != null)
        {
            HomologousRecombination.HrStatus hrStatus = HomologousRecombination.HrStatus.valueOf(chord.hrStatus().name());
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
                            .genes(GeneListUtil.genes(purple.reportableSomaticVariants(),
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

    private static FindingItem<MicrosatelliteStability> createMicrosatelliteStability(PurpleRecord purple,
            LinxRecord linx, DriverFindingList<GainDeletion> gainDeletions)
    {
        MicrosatelliteStability.MicrosatelliteStatus microsatelliteStatus =
                MicrosatelliteStability.MicrosatelliteStatus.valueOf(purple.characteristics()
                        .microsatelliteStatus()
                        .name());
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
                        .genes(GeneListUtil.genes(purple.reportableSomaticVariants(),
                                purple.reportableSomaticGainsDels(),
                                linx.germlineHomozygousDisruptions(),
                                Genes.MSI_GENES).stream().toList())
                        .build())
                .build();
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

    private static FindingsStatus purpleFindingsStatus(PurpleRecord purpleRecord)
    {
        return purpleRecord.fit().qc().status().equals(Set.of(PurpleQCStatus.PASS)) ? FindingsStatus.OK : FindingsStatus.NOT_AVAILABLE;
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
                .domainsKept(fusion.domainsKept())
                .domainsLost(fusion.domainsLost())
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

    private static FindingList<PharmocoGenotype> createPharmcoGenotypesFindings(@Nullable Set<PeachGenotype> peachGenotypes,
            boolean hasContamination)
    {
        if(peachGenotypes != null)
        {
            return FindingListBuilder.<PharmocoGenotype>builder()
                    .status(hasContamination ? FindingsStatus.NOT_RELIABLE : FindingsStatus.OK)
                    .findings(peachGenotypes.stream().map(o ->
                                    PharmocoGenotypeBuilder.builder()
                                            .findingKey(FindingKeys.pharmacoGenotype(o.gene(), o.allele()))
                                            .gene(o.gene())
                                            .allele(o.allele())
                                            .alleleCount(o.alleleCount())
                                            .function(o.function())
                                            .haplotype(o.haplotype())
                                            .linkedDrugs(o.linkedDrugs())
                                            .urlPrescriptionInfo(o.urlPrescriptionInfo())
                                            .build())
                            .sorted(PharmocoGenotype.COMPARATOR)
                            .toList())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingList();
        }
    }

    private static Map<String, DriverGene> driverGenesMap(@Nullable Path driverGeneTsv) throws IOException
    {
        return driverGeneTsv != null ? DriverGeneFile.read(driverGeneTsv)
                .stream()
                .collect(Collectors.toMap(DriverGene::gene, Function.identity())) : Map.of();
    }
}
