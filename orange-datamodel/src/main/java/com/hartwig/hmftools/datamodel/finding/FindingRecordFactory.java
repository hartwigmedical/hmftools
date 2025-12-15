package com.hartwig.hmftools.datamodel.finding;

import static com.hartwig.hmftools.datamodel.finding.DisruptionFactory.createGermlineDisruptions;
import static com.hartwig.hmftools.datamodel.finding.DisruptionFactory.createSomaticDisruptions;
import static com.hartwig.hmftools.datamodel.finding.GainDeletionFactory.germlineDriverGainDels;
import static com.hartwig.hmftools.datamodel.finding.GainDeletionFactory.somaticDriverGainDels;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Objects;

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
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory {
    public static final Logger LOGGER = LogManager.getLogger(FindingRecordFactory.class);

    @NotNull
    public static FindingRecord fromOrangeJsonWithTranscriptFile(@NotNull Path orangeJson, @NotNull Path clinicalTranscriptsTsv) throws IOException {
        try (Reader reader = Files.newBufferedReader(orangeJson)) {
            OrangeRecord orangeRecord = com.hartwig.hmftools.datamodel.OrangeJson.getInstance().read(reader);
            return fromOrangeRecordWithTranscriptFile(orangeRecord, clinicalTranscriptsTsv);
        }
    }

    @NotNull
    public static FindingRecord fromOrangeRecordWithTranscriptFile(@NotNull OrangeRecord orangeRecord, @NotNull Path clinicalTranscriptsTsv) throws IOException {
        ClinicalTranscriptsModel clinicalTranscriptsModel =
                ClinicalTranscriptFile.buildFromTsv(orangeRecord.refGenomeVersion(), clinicalTranscriptsTsv);
        return fromOrangeRecord(orangeRecord, clinicalTranscriptsModel);
    }

    @NotNull
    public static FindingRecord fromOrangeRecord(@NotNull OrangeRecord orangeRecord, @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel)
    {
        ImmutableFindingRecord.Builder builder = ImmutableFindingRecord.builder()
                .refGenomeVersion(orangeRecord.refGenomeVersion())
                .experimentType(orangeRecord.experimentType())
                .pipelineVersion(orangeRecord.pipelineVersion())
                .purpleFit(orangeRecord.purple().fit());

        builder = addPurpleFindings(builder, orangeRecord, clinicalTranscriptsModel);

        LinxRecord linx = orangeRecord.linx();
        boolean hasReliablePurity = orangeRecord.purple().fit().containsTumorCells();

        builder.driverSomaticDisruptions(createSomaticDisruptions(
                linx.reportableSomaticBreakends(),
                linx.allSomaticStructuralVariants(),
                linx.somaticDrivers(),
                hasReliablePurity));

        if(linx.allGermlineStructuralVariants() != null)
        {
            builder.driverGermlineDisruptions(
                    createGermlineDisruptions(
                            Objects.requireNonNull(linx.reportableGermlineBreakends()),
                            Objects.requireNonNull(linx.allGermlineStructuralVariants()),
                            Objects.requireNonNull(linx.germlineHomozygousDisruptions())));
        }

        CuppaData cuppa = orangeRecord.cuppa();

        if(cuppa != null)
        {
            builder.predictedTumorOrigin(ImmutablePredictedTumorOrigin.builder()
                    .findingKey(FindingKeys.predictedTumorOrigin(cuppa.bestPrediction().cancerType()))
                    .cancerType(cuppa.bestPrediction().cancerType())
                    .likelihood(cuppa.bestPrediction().likelihood())
                    .build());
        }

        ChordRecord chord = orangeRecord.chord();

        if(chord != null)
        {
            builder.homologousRecombination(ImmutableHomologousRecombination.builder()
                    .findingKey(FindingKeys.homologousRecombination(chord.hrStatus()))
                    .brca1Value(chord.brca1Value())
                    .brca2Value(chord.brca2Value())
                    .hrdValue(chord.hrdValue())
                    .hrStatus(chord.hrStatus())
                    .hrdType(chord.hrdType())
                    .build());
        }

        VirusInterpreterData virusInterpreter = orangeRecord.virusInterpreter();

        if(virusInterpreter != null)
        {
            builder.driverViruses(convertViruses(virusInterpreter.allViruses()));
        }

        builder.hlaFindings(HlaAlleleFactory.createHlaAllelesFinding(orangeRecord, hasReliablePurity));

        return builder.build();
    }

    @NotNull
    private static ImmutableFindingRecord.Builder addPurpleFindings(ImmutableFindingRecord.Builder builder, final OrangeRecord orangeRecord,
            final @Nullable ClinicalTranscriptsModel clinicalTranscriptsModel) {
        PurpleRecord purple = orangeRecord.purple();

        builder.driverSomaticSmallVariants(SmallVariantFactory.create(
                        DriverSource.SOMATIC, purple.reportableSomaticVariants(), orangeRecord.purple().somaticDrivers(),
                        clinicalTranscriptsModel))
                .driverSomaticGainDeletions(somaticDriverGainDels(purple.reportableSomaticGainsDels(), purple.somaticDrivers()))
                .driverSomaticFusions(orangeRecord.linx().reportableSomaticFusions().stream()
                        .map(o -> convertFusion(o, DriverSource.SOMATIC)).toList())
                .microsatelliteStability(
                        ImmutableMicrosatelliteStability.builder()
                                .findingKey(FindingKeys.microsatelliteStability(purple.characteristics().microsatelliteStatus()))
                                .microsatelliteStatus(purple.characteristics().microsatelliteStatus())
                                .microsatelliteIndelsPerMb(purple.characteristics().microsatelliteIndelsPerMb())
                                .build())
                .tumorMutationStatus(ImmutableTumorMutationStatus.builder()
                        .findingKey(FindingKeys.tumorMutationStatus(purple.characteristics().tumorMutationalBurdenStatus(),
                                purple.characteristics().tumorMutationalLoadStatus()))
                        .tumorMutationalBurdenPerMb(purple.characteristics().tumorMutationalBurdenPerMb())
                        .tumorMutationalBurdenStatus(purple.characteristics().tumorMutationalBurdenStatus())
                        .tumorMutationalLoad(purple.characteristics().tumorMutationalLoad())
                        .tumorMutationalLoadStatus(purple.characteristics().tumorMutationalLoadStatus())
                        .svTumorMutationalBurden(purple.characteristics().svTumorMutationalBurden())
                        .build());


        List<PurpleGainDeletion> reportableGermlineFullDels = orangeRecord.purple().reportableGermlineFullDels();
        List<PurpleLossOfHeterozygosity> reportableGermlineLohs = orangeRecord.purple().reportableGermlineLossOfHeterozygosities();
        List<PurpleDriver> purpleGermlineDrivers = orangeRecord.purple().germlineDrivers();

        if(reportableGermlineFullDels != null && reportableGermlineLohs != null && purpleGermlineDrivers != null) {
            builder.driverGermlineGainDeletions(germlineDriverGainDels(reportableGermlineFullDels, reportableGermlineLohs, purpleGermlineDrivers));
        }
        return builder;
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
}
