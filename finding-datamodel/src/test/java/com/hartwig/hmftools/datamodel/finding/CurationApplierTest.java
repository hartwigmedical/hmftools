package com.hartwig.hmftools.datamodel.finding;

import static com.hartwig.hmftools.datamodel.finding.TestFindingRecordFactory.createMinimalTestFindingRecordBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;

import org.junit.Test;

public class CurationApplierTest
{
    @Test
    public void applyWithNoCurationsReturnsUnchangedRecord()
    {
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-1", ReportedStatus.CANDIDATE))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .build();

        CurationRecord curation = new CurationRecord(List.of());
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.CANDIDATE, result.somaticSmallVariants().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationChangesSmallVariantToReported()
    {
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-1", ReportedStatus.CANDIDATE))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .build();

        DriverCuration dc = new DriverCuration("variant-1", 0L, "user1", true, "promote to reported");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.somaticSmallVariants().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationChangesSmallVariantToCandidate()
    {
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-1", ReportedStatus.REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .build();

        DriverCuration dc = new DriverCuration("variant-1", 0L, "user1", false, "demote to candidate");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.CANDIDATE, result.somaticSmallVariants().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationDoesNotAffectUnmatchedFindings()
    {
        SmallVariant variant1 = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-1", ReportedStatus.CANDIDATE))
                .build();
        SmallVariant variant2 = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-2", ReportedStatus.NOT_REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant1, variant2)))
                .build();

        DriverCuration dc = new DriverCuration("variant-1", 0L, "user1", true, "promote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.somaticSmallVariants().findings().get(0).reportedStatus());
        assertEquals(ReportedStatus.NOT_REPORTED, result.somaticSmallVariants().findings().get(1).reportedStatus());
    }

    @Test
    public void applyCurationToGainDeletion()
    {
        GainDeletion gd = TestFindingFactory.gainDeletionBuilder()
                .driver(driverFields("gd-1", ReportedStatus.NOT_REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticGainDeletions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(gd)))
                .build();

        DriverCuration dc = new DriverCuration("gd-1", 0L, "user1", true, "promote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.somaticGainDeletions().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationToDisruption()
    {
        Disruption disruption = TestFindingFactory.disruptionBuilder()
                .driver(driverFields("dis-1", ReportedStatus.REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticDisruptions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(disruption)))
                .build();

        DriverCuration dc = new DriverCuration("dis-1", 0L, "user1", false, "demote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.CANDIDATE, result.somaticDisruptions().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationToFusion()
    {
        Fusion fusion = TestFindingFactory.fusionBuilder()
                .driver(driverFields("fusion-1", ReportedStatus.CANDIDATE))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(fusion)))
                .build();

        DriverCuration dc = new DriverCuration("fusion-1", 0L, "user1", true, "promote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.fusions().findings().get(0).reportedStatus());
    }

    @Test
    public void applyCurationToGermlineList()
    {
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(driverFields("gl-variant-1", ReportedStatus.NOT_REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .germlineSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .build();

        DriverCuration dc = new DriverCuration("gl-variant-1", 0L, "user1", true, "promote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.germlineSmallVariants().findings().get(0).reportedStatus());
    }

    @Test
    public void applyPreservesOtherFieldsOnDriver()
    {
        DriverFields originalDriver = DriverFieldsBuilder.builder()
                .findingKey("variant-1")
                .driverSource(DriverSource.GERMLINE)
                .reportedStatus(ReportedStatus.CANDIDATE)
                .driverInterpretation(DriverInterpretation.LOW)
                .driverLikelihood(0.75)
                .build();
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(originalDriver)
                .gene("BRCA1")
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .build();

        DriverCuration dc = new DriverCuration("variant-1", 0L, "user1", true, "promote");
        CurationRecord curation = new CurationRecord(List.of(dc));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        SmallVariant resultVariant = result.somaticSmallVariants().findings().get(0);
        assertEquals(ReportedStatus.REPORTED, resultVariant.reportedStatus());
        assertEquals(DriverSource.GERMLINE, resultVariant.driverSource());
        assertEquals(DriverInterpretation.LOW, resultVariant.driverInterpretation());
        assertEquals(0.75, resultVariant.driverLikelihood(), 0.001);
        assertEquals("BRCA1", resultVariant.gene());
    }

    @Test
    public void applyMultipleCurationsAcrossDifferentLists()
    {
        SmallVariant variant = TestFindingFactory.variantBuilder()
                .driver(driverFields("variant-1", ReportedStatus.CANDIDATE))
                .build();
        Fusion fusion = TestFindingFactory.fusionBuilder()
                .driver(driverFields("fusion-1", ReportedStatus.REPORTED))
                .build();
        FindingRecord record = createMinimalTestFindingRecordBuilder()
                .somaticSmallVariants(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(variant)))
                .fusions(TestFindingFactory.buildDriverFindingsList(FindingsStatus.OK, List.of(fusion)))
                .build();

        CurationRecord curation = new CurationRecord(List.of(
                new DriverCuration("variant-1", 0L, "user1", true, "promote variant"),
                new DriverCuration("fusion-1", 0L, "user1", false, "demote fusion")));
        FindingRecord result = CurationApplier.applyCurations(record, curation);

        assertEquals(ReportedStatus.REPORTED, result.somaticSmallVariants().findings().get(0).reportedStatus());
        assertEquals(ReportedStatus.CANDIDATE, result.fusions().findings().get(0).reportedStatus());
    }

    private static DriverFields driverFields(String findingKey, ReportedStatus reportedStatus)
    {
        return DriverFieldsBuilder.builder()
                .findingKey(findingKey)
                .driverSource(DriverSource.SOMATIC)
                .reportedStatus(reportedStatus)
                .driverInterpretation(DriverInterpretation.HIGH)
                .build();
    }
}
