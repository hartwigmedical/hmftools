package com.hartwig.hmftools.finding.datamodel;

import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import org.jetbrains.annotations.NotNull;

@SuppressWarnings("unused")
public class TestFindingRecordFactory
{
    @NotNull
    public static FindingRecord createMinimalTestFindingRecord()
    {
        return createMinimalTestFindingRecordBuilder().build();
    }

    @NotNull
    public static FindingRecordBuilder createMinimalTestFindingRecordBuilder()
    {
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
                .predictedTumorOrigin(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.predictedTumorOriginBuilder()
                        .build()))
                .homologousRecombination(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.homologousRecombinationBuilder()
                        .build()))
                .microsatelliteStability(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.microsatelliteStabilityBuilder()
                        .build()))
                .tumorMutationalLoad(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.tumorMutationalLoadBuilder()
                        .build()))
                .tumorMutationalBurden(TestFindingFactory.buildFindingItem(FindingStatus.Status.OK, TestFindingFactory.tumorMutationalBurdenBuilder()
                        .build()))
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .germlineSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .somaticGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .germlineGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .somaticDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .germlineDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .chromosomeArmCopyNumbers(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of()))
                .viruses(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of()))
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of()))
                .pharmacoGenotypes(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of()));
    }
}
