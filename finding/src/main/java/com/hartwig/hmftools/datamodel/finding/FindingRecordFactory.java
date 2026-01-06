package com.hartwig.hmftools.datamodel.finding;

import static com.hartwig.hmftools.datamodel.finding.DisruptionFactory.createDisruptionsFindings;

import java.io.IOException;
import java.io.Reader;
import java.math.BigDecimal;
import java.math.RoundingMode;
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
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.clinicaltranscript.ClinicalTranscriptFile;
import com.hartwig.hmftools.datamodel.finding.clinicaltranscript.ClinicalTranscriptsModel;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.Genes;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory {

    @NotNull
    public static FindingRecord fromOrangeJsonWithTranscriptFile(@NotNull Path orangeJson, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException {
        try (Reader reader = Files.newBufferedReader(orangeJson)) {
            OrangeRecord orangeRecord = com.hartwig.hmftools.datamodel.OrangeJson.getInstance().read(reader);
            return fromOrangeRecord(orangeRecord, clinicalTranscriptsTsv, driverGeneTsv);
        }
    }

    @NotNull
    public static FindingRecord fromOrangeRecord(@NotNull OrangeRecord orangeRecord, @Nullable Path clinicalTranscriptsTsv,
            @Nullable Path driverGeneTsv) throws IOException {
        ClinicalTranscriptsModel clinicalTranscriptsModel = clinicalTranscriptsTsv != null ?
                ClinicalTranscriptFile.buildFromTsv(orangeRecord.refGenomeVersion(), clinicalTranscriptsTsv) : null;
        Map<String, DriverGene> driverGenes = driverGenesMap(driverGeneTsv);

        LinxRecord linx = orangeRecord.linx();
        boolean hasReliablePurity = orangeRecord.purple().fit().containsTumorCells();

        ImmutableFindingRecord.Builder builder = ImmutableFindingRecord.builder()
                .refGenomeVersion(orangeRecord.refGenomeVersion())
                .experimentType(orangeRecord.experimentType())
                .pipelineVersion(orangeRecord.pipelineVersion())
                .purpleFit(orangeRecord.purple().fit())
                .disruptions(createDisruptionsFindings(linx, hasReliablePurity))
                .fusions(createFusionsFindings(orangeRecord.linx()));

        addPurpleFindings(builder, orangeRecord, clinicalTranscriptsModel, driverGenes);

        return builder.predictedTumorOrigin(createPredictedTumorOrigin(orangeRecord.cuppa()))
                .homologousRecombination(createHomologousRecombination(orangeRecord.chord(), orangeRecord.purple(), linx))
                .viruses(createVirusFindings(orangeRecord.virusInterpreter()))
                .hla(HlaAlleleFactory.createHlaAllelesFindings(orangeRecord, hasReliablePurity))
                .pharmocoGenotypes(createPharmcoGenotypesFindings(orangeRecord.peach()))
                .build();
    }

    private static void addPurpleFindings(ImmutableFindingRecord.Builder builder, final OrangeRecord orangeRecord,
            final @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel, @NotNull Map<String, DriverGene> driverGenes) {
        PurpleRecord purple = orangeRecord.purple();

        FindingsStatus findingsStatus = purpleFindingsStatus(purple);

        builder.smallVariants(SmallVariantFactory.smallVariantFindings(purple, findingsStatus, clinicalTranscriptsModel, driverGenes))
                .gainDeletions(GainDeletionFactory.gainDeletionFindings(purple, findingsStatus))
                .microsatelliteStability(createMicrosatelliteStability(purple, orangeRecord.linx()))
                .tumorMutationStatus(createTumorMutationStatus(purple));

        builder.chromosomeArmCopyNumbers(
                ChromosomeArmCopyNumberFactory.extractCnPerChromosomeArm(purple.allSomaticCopyNumbers(), orangeRecord.refGenomeVersion()));
    }

    @NotNull
    private static CharacteristicsFinding<TumorMutationStatus> createTumorMutationStatus(@NotNull PurpleRecord purple) {
        return ImmutableCharacteristicsFinding.<TumorMutationStatus>builder()
                .status(FindingsStatus.OK)
                .finding(ImmutableTumorMutationStatus.builder()
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

    @NotNull
    private static CharacteristicsFinding<PredictedTumorOrigin> createPredictedTumorOrigin(CuppaData cuppa) {
        if (cuppa != null) {
            return ImmutableCharacteristicsFinding.<PredictedTumorOrigin>builder()
                    .status(FindingsStatus.OK)
                    .finding( ImmutablePredictedTumorOrigin.builder()
                            .findingKey(FindingKeys.predictedTumorOrigin(cuppa.bestPrediction().cancerType()))
                            .cancerType(cuppa.bestPrediction().cancerType())
                            .likelihood(cuppa.bestPrediction().likelihood())
                            .build())
                    .build();
        } else {
            return ImmutableCharacteristicsFinding.<PredictedTumorOrigin>builder()
                    .status(FindingsStatus.NOT_AVAILABLE)
                    .build();
        }
    }

    @NotNull
    private static CharacteristicsFinding<HomologousRecombination> createHomologousRecombination(@Nullable ChordRecord chord,
            @NotNull PurpleRecord purple,
            @NotNull LinxRecord linx) {
        if (chord != null) {
            return ImmutableCharacteristicsFinding.<HomologousRecombination>builder()
                    .status(FindingsStatus.OK)
                    .finding(ImmutableHomologousRecombination.builder()
                            .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                            .brca1Value(chord.brca1Value())
                            .brca2Value(chord.brca2Value())
                            .hrdValue(chord.hrdValue())
                            .hrStatus(chord.hrStatus())
                            .hrdType(chord.hrdType())
                            .lohCopyNumbers(createGeneCopyNumbers(purple, Genes.HRD_GENES ))
                            .genes(GeneListUtil.genes(purple.reportableSomaticVariants(),
                                    purple.reportableSomaticGainsDels(),
                                    linx.germlineHomozygousDisruptions(),
                                    Genes.HRD_GENES))
                            .build())
                    .build();
        } else {
            return ImmutableCharacteristicsFinding.<HomologousRecombination>builder()
                    .status(FindingsStatus.NOT_AVAILABLE)
                    .build();
        }
    }

    @NotNull
    private static CharacteristicsFinding<MicrosatelliteStability> createMicrosatelliteStability(@NotNull PurpleRecord purple, @NotNull LinxRecord linx) {
        return ImmutableCharacteristicsFinding.<MicrosatelliteStability>builder()
                .status(FindingsStatus.OK)
                .finding(ImmutableMicrosatelliteStability.builder()
                        .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                        .microsatelliteStatus(purple.characteristics().microsatelliteStatus())
                        .microsatelliteIndelsPerMb(purple.characteristics().microsatelliteIndelsPerMb())
                        .lohCopyNumbers(createGeneCopyNumbers(purple, Genes.MSI_GENES))
                        .genes(GeneListUtil.genes(purple.reportableSomaticVariants(),
                                purple.reportableSomaticGainsDels(),
                                linx.germlineHomozygousDisruptions(),
                                Genes.MSI_GENES))
                        .build())
                .build();
    }

    @NotNull
    private static List<LOHCopyNumbers> createGeneCopyNumbers(@NotNull PurpleRecord purpleRecord, Set<String> geneNames) {
        List<PurpleGeneCopyNumber> suspectGeneCopyNumbersWithLOH = purpleRecord.suspectGeneCopyNumbersWithLOH();
        boolean hasReliablePurity = purpleRecord.fit().containsTumorCells();
        return suspectGeneCopyNumbersWithLOH.stream()
                .filter(x -> geneNames.contains(x.gene()))
                .map(lohGene -> ImmutableLOHCopyNumbers.builder()
                        .findingKey(FindingKeys.lohCopyNumber(lohGene))
                        .location(lohGene.chromosome() + lohGene.chromosomeBand())
                        .gene(lohGene.gene())
                        .tumorCopies(hasReliablePurity ? toInteger(lohGene.minCopyNumber()) : null)
                        .tumorMinorAlleleCopies(hasReliablePurity ? toInteger(lohGene.minMinorAlleleCopyNumber()) : null)
                        .build())
                .sorted(Comparator.comparing(LOHCopyNumbers::gene))
                .collect(Collectors.toList());
    }

    private static Integer toInteger(@Nullable Double value) {
        return value != null ? BigDecimal.valueOf(value).setScale(0, RoundingMode.HALF_EVEN) // match DecimalFormat
                .intValueExact() : null;
    }

    private static FindingsStatus purpleFindingsStatus(PurpleRecord purpleRecord) {
        return purpleRecord.fit().qc().status().equals(Set.of(PurpleQCStatus.PASS)) ? FindingsStatus.OK : FindingsStatus.NOT_AVAILABLE;
    }

    public static DriverFindings<Fusion> createFusionsFindings(@NotNull LinxRecord linx) {
        return ImmutableDriverFindings.<Fusion>builder()
                .status(FindingsStatus.OK)
                .all(linx.reportableSomaticFusions().stream()
                .map(o -> convertFusion(o, DriverSource.SOMATIC)).toList())
                .build();
    }

    public static Fusion convertFusion(LinxFusion fusion, DriverSource sampleType)
    {
        DriverInterpretation driverInterpretation = toDriverInterpretation(fusion.driverLikelihood());

        return ImmutableFusion.builder()
                .findingKey(FindingKeys.fusion(sampleType, fusion))
                .driverSource(sampleType)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(driverInterpretation)
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

    @NotNull
    private static DriverInterpretation toDriverInterpretation(@NotNull FusionLikelihoodType likelihood)
    {
        return switch (likelihood) {
            case HIGH -> DriverInterpretation.HIGH;
            case LOW -> DriverInterpretation.LOW;
            case NA -> DriverInterpretation.UNKNOWN;
        };
    }

    @NotNull
    private static DriverFindings<Virus> createVirusFindings(@Nullable VirusInterpreterData virusInterpreter){
        if (virusInterpreter != null) {
            return ImmutableDriverFindings.<Virus>builder()
                    .status(FindingsStatus.OK)
                    .all(convertViruses(virusInterpreter.allViruses()))
                    .build();

        } else {
            return ImmutableDriverFindings.<Virus>builder()
                    .status(FindingsStatus.NOT_APPLICABLE)
                    .build();
        }
    }

    @NotNull
    private static List<? extends Virus> convertViruses(List<VirusInterpreterEntry> viruses)
    {
        return viruses.stream()
                .map(v -> ImmutableVirus.builder()
                        .findingKey(FindingKeys.virus(v))
                        .driverSource(DriverSource.SOMATIC)
                        .reportedStatus(v.reported() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED)
                        .driverInterpretation(virusDriverInterpretation(v.driverLikelihood()))
                        .interpreterEntry(v)
                        .build()
                ).toList();
    }

    @NotNull
    private static DriverInterpretation virusDriverInterpretation(@NotNull VirusLikelihoodType virusLikelihoodType)
    {
        return switch(virusLikelihoodType)
        {
            case LOW -> DriverInterpretation.LOW;
            case HIGH -> DriverInterpretation.HIGH;
            case UNKNOWN -> DriverInterpretation.UNKNOWN;
        };
    }

    private static Findings<PharmocoGenotype> createPharmcoGenotypesFindings(@Nullable Set<PeachGenotype> peachGenotypes) {
        if(peachGenotypes != null)
        {
            return ImmutableFindings.<PharmocoGenotype>builder()
                    .status(FindingsStatus.OK)
                    .all(peachGenotypes.stream().map(o ->
                            ImmutablePharmocoGenotype.builder()
                                    .findingKey(FindingKeys.pharmacoGenotype(o.gene(), o.allele()))
                                    .gene(o.gene())
                                    .allele(o.allele())
                                    .alleleCount(o.alleleCount())
                                    .function(o.function())
                                    .haplotype(o.haplotype())
                                    .linkedDrugs(o.linkedDrugs())
                                    .urlPrescriptionInfo(o.urlPrescriptionInfo())
                                    .build()).toList())
                    .build();
        }
        else
        {
            return ImmutableFindings.<PharmocoGenotype>builder()
                    .status(FindingsStatus.NOT_AVAILABLE)
                    .build();
        }
    }

    private static Map<String, DriverGene> driverGenesMap(@Nullable Path driverGeneTsv) throws IOException {
        return driverGeneTsv != null ? DriverGeneFile.read(driverGeneTsv).stream().collect(Collectors.toMap(DriverGene::gene, Function.identity())) : Map.of();
    }
}
