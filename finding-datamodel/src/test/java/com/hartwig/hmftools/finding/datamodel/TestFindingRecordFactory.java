package com.hartwig.hmftools.finding.datamodel;

import java.time.LocalDate;
import java.util.List;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("unused")
public class TestFindingRecordFactory {
    @NotNull
    public static FindingRecord createMinimalTestFindingRecord() {
        return createMinimalTestFindingRecordBuilder().build();
    }

    @NotNull
    public static FindingRecordBuilder createMinimalTestFindingRecordBuilder() {
        return FindingRecordBuilder.builder()
                .version("")
                .metaProperties(MetaPropertiesBuilder.builder()
                        .sampleId("")
                        .samplingDate(LocalDate.of(2026, 1, 1))
                        .sequencingScope(SequencingScope.WHOLE_GENOME)
                        .refGenomeVersion(RefGenomeVersion.V37)
                        .build())
                .qc(TestFindingFactory.qcBuilder().build())
                .purityPloidyFit(TestFindingFactory.purityPloidyFitBuilder().build())
                .predictedTumorOrigins(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of(TestFindingFactory.predictedTumorOriginBuilder().build())))
                .homologousRecombination(TestFindingFactory.buildFindingItem(FindingsStatus.OK, TestFindingFactory.homologousRecombinationBuilder().build()))
                .microsatelliteStability(TestFindingFactory.buildFindingItem(FindingsStatus.OK, TestFindingFactory.microsatelliteStabilityBuilder().build()))
                .tumorMutationStatus(TestFindingFactory.buildFindingItem(FindingsStatus.OK, TestFindingFactory.mutationStatusBuilder().build()))
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .germlineSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .somaticGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .germlineGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .somaticDisruptions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .germlineDisruptions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .viruses(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of()))
                .chromosomeArmCopyNumbers(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .pharmacoGenotypes(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .visualisationFiles(VisualisationFilesBuilder.builder()
                        .purpleInputPlot("")
                        .purpleFinalCircosPlot("")
                        .purpleClonalityPlot("")
                        .purpleCopyNumberPlot("")
                        .purpleVariantCopyNumberPlot("")
                        .purplePurityRangePlot("")
                        .purpleKataegisPlot("")
                        .qseePlot("")
                        .linxDriverPlots(List.of())
                        .sageVisualisations(List.of())
                        .build());
    }
}
