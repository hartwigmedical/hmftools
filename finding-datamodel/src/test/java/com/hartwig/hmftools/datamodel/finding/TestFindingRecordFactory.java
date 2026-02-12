package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

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
                        .experimentType(ExperimentType.WHOLE_GENOME)
                        .refGenomeVersion(OrangeRefGenomeVersion.V37)
                        .build())
                .purityPloidyFit(TestFindingFactory.purityPloidyFitBuilder().build())
                .predictedTumorOrigin(TestFindingFactory.buildFindingItem(FindingsStatus.OK, TestFindingFactory.predictedTumorOriginBuilder().build()))
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
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .pharmocoGenotypes(TestFindingFactory.buildFindingsList(FindingsStatus.OK, List.of()))
                .visualisationFiles(VisualisationFilesBuilder.builder()
                        .tumorBqrPlot("")
                        .purpleInputPlot("")
                        .purpleFinalCircosPlot("")
                        .purpleClonalityPlot("")
                        .purpleCopyNumberPlot("")
                        .purpleVariantCopyNumberPlot("")
                        .purplePurityRangePlot("")
                        .purpleKataegisPlot("")
                        .qseePlot("")
                        .linxDriverPlots(List.of())
                        .sageVisualisations(Map.of())
                        .build());
    }
}
