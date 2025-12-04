package com.hartwig.hmftools.datamodel.finding;

import static com.hartwig.hmftools.datamodel.finding.DisruptionFactory.createDisruptions;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// to reduce duplication, the findings are collected from
// various part of the orange record
public class FindingRecordFactory
{
    public static final Logger LOGGER = LogManager.getLogger(FindingRecordFactory.class);

    public static FindingRecord fromOrangeRecord(OrangeRecord orangeRecord)
    {

        PurpleRecord purple = orangeRecord.purple();

        ImmutableFindingRecord.Builder builder = ImmutableFindingRecord.builder()
                .driverSomaticSmallVariants(SmallVariantFactory.create(
                        FindingKeys.SampleType.SOMATIC, purple.reportableSomaticVariants(), orangeRecord.purple().somaticDrivers()))
                .driverSomaticFusions(orangeRecord.linx().reportableSomaticFusions().stream()
                        .map(o -> convertFusion(o, FindingKeys.SampleType.SOMATIC)).toList())
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

        List<PurpleVariant> germlineVariants = orangeRecord.purple().reportableGermlineVariants();
        List<PurpleDriver> germlineDrivers = orangeRecord.purple().germlineDrivers();
        if (germlineVariants != null && germlineDrivers != null) {
            builder.driverGermlineSmallVariants(SmallVariantFactory.create(
                    FindingKeys.SampleType.GERMLINE, germlineVariants, germlineDrivers));
        }

        LinxRecord linx = orangeRecord.linx();
        boolean hasReliablePurity = orangeRecord.purple().fit().containsTumorCells();

        builder.driverSomaticDisruptions(createDisruptions(FindingKeys.SampleType.SOMATIC, linx.reportableSomaticBreakends(),
                linx.allSomaticStructuralVariants(), hasReliablePurity));

        if(linx.allGermlineStructuralVariants() != null)
        {
            builder.driverGermlineDisruptions(createDisruptions(FindingKeys.SampleType.GERMLINE,
                    Objects.requireNonNull(linx.reportableGermlineBreakends()),
                    Objects.requireNonNull(linx.allGermlineStructuralVariants()), hasReliablePurity));
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
            builder.driverViruses(virusInterpreter.reportableViruses().stream()
                    .map(v -> ImmutableVirus.builder()
                            .findingKey(FindingKeys.virus(v))
                            .reportedStatus(ReportedStatus.REPORTED)
                            .driverInterpretation(virusDriverInterpretation(v.driverLikelihood()))
                            .interpreterEntry(v)
                            .build()
                    ).toList());
        }

        return builder.build();
    }

    public static Fusion convertFusion(LinxFusion fusion, FindingKeys.SampleType sampleType)
    {
        DriverInterpretation driverInterpretation = toDriverInterpretation(fusion.driverLikelihood());

        return ImmutableFusion.builder()
                .findingKey(FindingKeys.fusion(sampleType, fusion))
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
            case LOW, NA -> DriverInterpretation.LOW;
        };
    }

    @NotNull
    private static DriverInterpretation virusDriverInterpretation(@NotNull VirusLikelihoodType virusLikelihoodType)
    {
        return switch(virusLikelihoodType)
        {
            case LOW -> DriverInterpretation.LOW;
            case HIGH -> DriverInterpretation.HIGH;
            default -> throw new IllegalStateException("Unexpected virus likelihood type: " + virusLikelihoodType);
        };
    }
}
