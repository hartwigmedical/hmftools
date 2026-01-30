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
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.finding.*;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.finding.clinicaltranscript.ClinicalTranscriptsModel;
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
        boolean hasContainmination = purple.fit().qc().status().contains(PurpleQCStatus.FAIL_CONTAMINATION);

        FindingRecordBuilder builder = FindingRecordBuilder.builder()
                .metaProperties(MetaPropertiesBuilder.builder()
                        .version("1.0")
                        .refGenomeVersion(orangeRecord.refGenomeVersion())
                        .experimentType(orangeRecord.experimentType())
                        .pipelineVersion(orangeRecord.pipelineVersion())
                        .build())
                .somaticDisruptions(createSomaticDisruptions(hasReliablePurity, linx))
                .germlineDisruptions(createGermlineDisruptions(orangeRecord.refSample() != null, linx))
                .fusions(createFusionsFindings(orangeRecord.linx()));

        DriverFindingList<GainDeletion> somaticGainDeletions = addPurpleFindings(builder, orangeRecord, clinicalTranscriptsModel, driverGenes);

        return builder.predictedTumorOrigin(createPredictedTumorOrigin(orangeRecord.cuppa()))
                .homologousRecombination(createHomologousRecombination(orangeRecord.chord(), purple, linx, somaticGainDeletions))
                .viruses(createVirusFindings(orangeRecord.virusInterpreter()))
                .hla(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, hasReliablePurity))
                .pharmocoGenotypes(createPharmcoGenotypesFindings(orangeRecord.peach(), hasContainmination))
                .build();
    }

    // return the gain deletions cause they are needed by HRD, will see if we can find a better way
    private static DriverFindingList<GainDeletion> addPurpleFindings(FindingRecordBuilder builder, final OrangeRecord orangeRecord,
            final @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel, Map<String, DriverGene> driverGenes)
    {
        boolean hasRefSample = orangeRecord.refSample() != null;

        PurpleRecord purple = orangeRecord.purple();

        FindingsStatus findingsStatus = purpleFindingsStatus(purple);

        PurpleFit purpleFit = purple.fit();

        builder.purityPloidyFit(PurityPloidyFitBuilder.builder()
                .qc(PurityPloidyFitQcBuilder.builder()
                        .status(purpleFit.qc().status())
                        .germlineAberrations(purpleFit.qc().germlineAberrations())
                        .amberMeanDepth(purpleFit.qc().amberMeanDepth())
                        .contamination(purpleFit.qc().contamination())
                        .totalCopyNumberSegments(purpleFit.qc().totalCopyNumberSegments())
                        .unsupportedCopyNumberSegments(purpleFit.qc().unsupportedCopyNumberSegments())
                        .deletedGenes(purpleFit.qc().deletedGenes())
                        .build())
                .fittedPurityMethod(purpleFit.fittedPurityMethod())
                .purity(purpleFit.purity())
                .minPurity(purpleFit.minPurity())
                .maxPurity(purpleFit.maxPurity())
                .ploidy(purpleFit.ploidy())
                .minPloidy(purpleFit.minPloidy())
                .maxPloidy(purpleFit.maxPloidy())
                .build());

        DriverFindingList<GainDeletion> somaticGainDeletions =
                GainDeletionFactory.somaticGainDeletionFindings(orangeRecord.refGenomeVersion(), findingsStatus, purple);

        builder.somaticSmallVariants(SmallVariantFactory.somaticSmallVariantFindings(purple, findingsStatus, clinicalTranscriptsModel, driverGenes))
                .germlineSmallVariants(SmallVariantFactory.germlineSmallVariantFindings(hasRefSample, purple, clinicalTranscriptsModel, driverGenes))
                .somaticGainDeletions(somaticGainDeletions)
                .germlineGainDeletions(GainDeletionFactory.germlineGainDeletionFindings(hasRefSample, orangeRecord.refGenomeVersion(), purple))
                .microsatelliteStability(createMicrosatelliteStability(purple, orangeRecord.linx(), somaticGainDeletions))
                .tumorMutationStatus(createTumorMutationStatus(purple));

        return somaticGainDeletions;
    }

    private static FindingItem<TumorMutationStatus> createTumorMutationStatus(PurpleRecord purple)
    {
        return FindingItemBuilder.<TumorMutationStatus>builder()
                .status(FindingsStatus.OK)
                .finding(TumorMutationStatusBuilder.builder()
                        .findingKey(FindingKeys.tumorMutationStatus(purple.characteristics().tumorMutationalBurdenStatus(),
                                purple.characteristics().tumorMutationalLoadStatus()))
                        .tumorMutationalBurdenPerMb(purple.characteristics().tumorMutationalBurdenPerMb())
                        .tumorMutationalBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus())
                        .tumorMutationalLoad(purple.characteristics().tumorMutationalLoad())
                        .tumorMutationalLoadStatus(purple.characteristics().tumorMutationalLoadStatus())
                        .svTumorMutationalBurden(purple.characteristics().svTumorMutationalBurden())
                        .build())
                .build();
    }

    private static FindingItem<PredictedTumorOrigin> createPredictedTumorOrigin(@Nullable CuppaData cuppa)
    {
        if(cuppa != null)
        {
            return FindingItemBuilder.<PredictedTumorOrigin>builder()
                    .status(FindingsStatus.OK)
                    .finding(PredictedTumorOriginBuilder.builder()
                            .findingKey(FindingKeys.predictedTumorOrigin(cuppa.bestPrediction().cancerType()))
                            .cancerType(cuppa.bestPrediction().cancerType())
                            .likelihood(cuppa.bestPrediction().likelihood())
                            .build())
                    .build();
        }
        else
        {
            return FindingUtil.notAvailableFindingItem();
        }
    }

    private static FindingItem<HomologousRecombination> createHomologousRecombination(@Nullable ChordRecord chord,
            PurpleRecord purple,
            LinxRecord linx,
            DriverFindingList<GainDeletion> gainDeletions)
    {
        if(chord != null)
        {
            return FindingItemBuilder.<HomologousRecombination>builder()
                    .status(FindingsStatus.OK)
                    .finding(HomologousRecombinationBuilder.builder()
                            .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                            .brca1Value(chord.brca1Value())
                            .brca2Value(chord.brca2Value())
                            .hrdValue(chord.hrdValue())
                            .hrStatus(chord.hrStatus())
                            .hrdType(chord.hrdType())
                            .lohCopyNumbers(filterLohGainDeletions(gainDeletions, Genes.HRD_GENES))
                            .genes(GeneListUtil.genes(purple.somaticVariants(),
                                    purple.somaticGainsDels(),
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
        return FindingItemBuilder.<MicrosatelliteStability>builder()
                .status(FindingsStatus.OK)
                .finding(MicrosatelliteStabilityBuilder.builder()
                        .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                        .microsatelliteStatus(purple.characteristics().microsatelliteStatus())
                        .microsatelliteIndelsPerMb(purple.characteristics().microsatelliteIndelsPerMb())
                        .lohCopyNumbers(filterLohGainDeletions(gainDeletions, Genes.MSI_GENES))
                        .genes(GeneListUtil.genes(purple.somaticVariants(),
                                purple.somaticGainsDels(),
                                linx.germlineHomozygousDisruptions(),
                                Genes.MSI_GENES).stream().toList())
                        .build())
                .build();
    }

    private static List<GainDeletion> filterLohGainDeletions(DriverFindingList<GainDeletion> gainDeletions, Set<String> geneNames)
    {
        return gainDeletions.findings().stream()
                .filter(x -> geneNames.contains(x.gene()))
                .filter(x -> x.type() == GainDeletion.Type.SOMATIC_LOH)
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
                .findings(linx.fusions().stream()
                        .map(o -> convertFusion(o, DriverSource.SOMATIC)).sorted(Fusion.COMPARATOR).toList())
                .build();
    }

    public static Fusion convertFusion(LinxFusion fusion, DriverSource sampleType)
    {
        DriverInterpretation driverInterpretation = toDriverInterpretation(fusion.driverLikelihood());

        boolean isDriverGene = !fusion.unreportedReasons().contains(LinxUnreportableReason.NOT_KNOWN);

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
                .reportedType(fusion.reportedType())
                .unreportedReasons(fusion.unreportedReasons())
                .phased(fusion.phased())
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
            case NA -> DriverInterpretation.LOW;
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
                        .qcStatus(v.qcStatus())
                        .integrations(v.integrations())
                        .interpretation(v.interpretation())
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
