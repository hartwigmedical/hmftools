package com.hartwig.hmftools.finding.datamodel;

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
                .metaProperties(MetaPropertiesBuilder.builder()
                        .sequencingScope(SequencingScope.WHOLE_GENOME)
                        .pipelineVersion("3.0")
                        .version("0.1.0")
                        .refGenomeVersion(RefGenomeVersion.V37)
                        .build())
                .purityPloidyFit(TestFindingFactory.purityPloidyFitBuilder().build())
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
                .predictedTumorOrigins(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .chromosomeArmCopyNumbers(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .pharmocoGenotypes(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
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
