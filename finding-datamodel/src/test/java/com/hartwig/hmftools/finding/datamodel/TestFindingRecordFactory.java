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
                .predictedTumorOrigin(TestFindingFactory.buildFindingItem(FindingStatus.ResultStatus.OK, TestFindingFactory.predictedTumorOriginBuilder()
                        .build()))
                .homologousRecombination(TestFindingFactory.buildFindingItem(FindingStatus.ResultStatus.OK, TestFindingFactory.homologousRecombinationBuilder()
                        .build()))
                .microsatelliteStability(TestFindingFactory.buildFindingItem(FindingStatus.ResultStatus.OK, TestFindingFactory.microsatelliteStabilityBuilder()
                        .build()))
                .tumorMutationalLoad(TestFindingFactory.buildFindingItem(FindingStatus.ResultStatus.OK, TestFindingFactory.tumorMutationalLoadBuilder()
                        .build()))
                .tumorMutationalBurden(TestFindingFactory.buildFindingItem(FindingStatus.ResultStatus.OK, TestFindingFactory.tumorMutationalBurdenBuilder()
                        .build()))
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .germlineSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .somaticGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .germlineGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .somaticDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .germlineDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .viruses(TestFindingFactory.buildDriverFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingStatus.ResultStatus.OK, List.of()))
                .pharmacoGenotypes(TestFindingFactory.buildFindingsList(FindingStatus.ResultStatus.OK, List.of()));
    }
}
