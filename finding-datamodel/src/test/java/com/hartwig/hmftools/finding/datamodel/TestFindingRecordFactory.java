package com.hartwig.hmftools.finding.datamodel;

import java.time.LocalDate;
import java.util.List;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.driver.DriverInterpretation;
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
                .metaProperties(metaPropertiesBuilder().build())
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
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.smallVariantBuilder()
                        .build())))
                .germlineSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.smallVariantBuilder()
                        .build())))
                .somaticGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.gainDeletionBuilder()
                        .build())))
                .germlineGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.gainDeletionBuilder()
                        .build())))
                .somaticDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.disruptionBuilder()
                        .build())))
                .germlineDisruptions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.disruptionBuilder()
                        .build())))
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.fusionBuilder()
                        .build())))
                .chromosomeArmCopyNumbers(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.chromosomeArmCopyNumberBuilder()
                        .build())))
                .viruses(TestFindingFactory.buildDriverFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.virusBuilder(true, DriverInterpretation.HIGH)
                        .build())))
                .hlaAlleles(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.hlaAlleleBuilder()
                        .build())))
                .pharmacoGenotypes(TestFindingFactory.buildFindingsList(FindingStatus.Status.OK, List.of(TestFindingFactory.pharmacoGenotypeBuilder()
                        .build())));
    }

    public static MetaPropertiesBuilder metaPropertiesBuilder()
    {
        return MetaPropertiesBuilder.builder()
                .sampleId("")
                .samplingDate(LocalDate.of(2026, 1, 1))
                .sequencingScope(SequencingScope.WHOLE_GENOME)
                .refGenomeVersion(RefGenomeVersion.V37)
                .potentialHRDGenes(new TreeSet<>())
                .potentialMSIGenes(new TreeSet<>());
    }
}
