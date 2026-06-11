package com.hartwig.hmftools.finding.util;

import static com.hartwig.hmftools.finding.datamodel.finding.FindingStatus.Issue.LOW_PURITY;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecordBuilder;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.Qc;
import com.hartwig.hmftools.finding.datamodel.SequencingScope;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.TestFindingRecordFactory;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.IFindingList;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class LowPurityTransformerTest
{
    private final static double EPSILON = 0.00001;

    @Test
    public void testLowPurityThresholdsWGS()
    {
        FindingRecord original = findingRecordBuilder(0.01, SequencingScope.WHOLE_GENOME).build();
        FindingRecord converted = LowPurityTransformer.transform(original);

        assertFindingListOk(converted.somaticSmallVariants());
        assertFindingListOk(converted.germlineSmallVariants());
        assertFindingListLowPurity(converted.somaticDisruptions(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.germlineDisruptions());
        assertFindingListLowPurity(converted.somaticGainDeletions(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.germlineGainDeletions());
        assertFindingListOk(converted.fusions());
        assertFindingListLowPurity(converted.viruses(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.chromosomeArmCopyNumbers());
        assertFindingListOk(converted.hlaAlleles(), LowPurityTransformer.MIN_PURITY_20_PCT, Set.of(LOW_PURITY));
        assertFindingListOk(converted.pharmacoGenotypes());

        assertFindingItemOk(converted.predictedTumorOrigin());
        assertFindingItemLowPurity(converted.homologousRecombination(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingItemLowPurity(converted.microsatelliteStability(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingItemLowPurity(converted.tumorMutationalBurden(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingItemLowPurity(converted.tumorMutationalLoad(), LowPurityTransformer.MIN_PURITY_20_PCT);
    }

    @Test
    public void testLowPurityThresholdsTargeted()
    {
        FindingRecord original = findingRecordBuilder(0.01, SequencingScope.TARGETED).build();
        FindingRecord converted = LowPurityTransformer.transform(original);

        assertFindingListLowPurity(converted.somaticSmallVariants(), LowPurityTransformer.MIN_PURITY_10_PCT);
        assertFindingListOk(converted.germlineSmallVariants());
        assertFindingListLowPurity(converted.somaticDisruptions(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.germlineDisruptions());
        assertFindingListLowPurity(converted.somaticGainDeletions(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.germlineGainDeletions());
        assertFindingListOk(converted.fusions());
        assertFindingListLowPurity(converted.viruses(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingListOk(converted.chromosomeArmCopyNumbers());
        assertFindingListOk(converted.hlaAlleles(), LowPurityTransformer.MIN_PURITY_20_PCT, Set.of(LOW_PURITY));
        assertFindingListOk(converted.pharmacoGenotypes());

        assertFindingItemOk(converted.predictedTumorOrigin());
        assertFindingItemLowPurity(converted.homologousRecombination(), LowPurityTransformer.MIN_PURITY_30_PCT);
        assertFindingItemLowPurity(converted.microsatelliteStability(), LowPurityTransformer.MIN_PURITY_20_PCT);
        assertFindingItemLowPurity(converted.tumorMutationalBurden(), LowPurityTransformer.MIN_PURITY_10_PCT);
        assertFindingItemLowPurity(converted.tumorMutationalLoad(), LowPurityTransformer.MIN_PURITY_10_PCT);
    }

    @Test
    public void testTransformHLALowPurity()
    {
        FindingRecord original = findingRecordBuilder(0.01, SequencingScope.WHOLE_GENOME)
                .hlaAlleles(FindingListBuilder.<HlaAllele>builder()
                        .status(TestFindingFactory.findingsStatus(FindingStatus.Status.OK))
                        .findings(List.of(TestFindingFactory.hlaAlleleBuilder().build()))
                        .build())
                .build();

        FindingRecord converted = LowPurityTransformer.transform(original);

        assertQc(converted.qc());
        assertHLA(converted.hlaAlleles());
    }

    private void assertHLA(FindingList<HlaAllele> findingList)
    {
        FindingStatus findingStatus = findingList.status();
        assertFindingStatus(findingStatus, FindingStatus.Status.OK, Set.of(), Set.of(LOW_PURITY));
        assertPurityThreshold(findingList, LowPurityTransformer.MIN_PURITY_20_PCT);
        List<HlaAllele> hlaAlleles = findingList.findings();
        assertFalse(hlaAlleles.isEmpty());
        for(HlaAllele hlaAllele : hlaAlleles)
        {
            assertNull(hlaAllele.tumorCopyNumber());
        }
    }

    @Test
    public void testTransformGainDeletionsNoLowPurity()
    {
        FindingRecord original = findingRecordBuilder(0.4, SequencingScope.WHOLE_GENOME)
                .somaticGainDeletions(TestFindingFactory.driverFindingsBuilder(List.of(TestFindingFactory.gainDeletionBuilder().build()))
                        .build())
                .build();
        FindingRecord converted = LowPurityTransformer.transform(original);

        DriverFindingList<?> gainDeletions = converted.somaticGainDeletions();
        assertFindingStatusOk(gainDeletions.status());
        assertEquals(1, gainDeletions.findings().size());
    }

    @Test
    public void testTransformGainDeletionsWgsLowPurity()
    {
        FindingRecord original = findingRecordBuilder(0.01, SequencingScope.WHOLE_GENOME)
                .somaticGainDeletions(TestFindingFactory.driverFindingsBuilder(List.of(TestFindingFactory.gainDeletionBuilder().build()))
                        .build())
                .build();
        FindingRecord converted = LowPurityTransformer.transform(original);

        DriverFindingList<?> gainDeletions = converted.somaticGainDeletions();
        assertFindingStatusLowPurity(gainDeletions.status());
        assertTrue(gainDeletions.findings().isEmpty());
    }

    @Test
    public void testTransformGainDeletionsTargetedLowPurityNoGainsNoDeletions()
    {
        testTransformGainDeletionsTargetedLowPurity(0.1, LowPurityTransformer.MIN_PURITY_20_PCT, FindingStatus.Status.NOT_RELIABLE, Set.of(LOW_PURITY), Set.of(), 0);
    }

    @Test
    public void testTransformGainDeletionsTargetedLowPurityNoDeletions()
    {
        testTransformGainDeletionsTargetedLowPurity(0.2, LowPurityTransformer.MIN_PURITY_30_PCT, FindingStatus.Status.OK, Set.of(), Set.of(LOW_PURITY), 1);
    }

    @Test
    public void testTransformGainDeletionsTargetedLowPurity()
    {
        testTransformGainDeletionsTargetedLowPurity(0.3, LowPurityTransformer.MIN_PURITY_20_PCT, FindingStatus.Status.OK, Set.of(), Set.of(), 2);
    }

    private void testTransformGainDeletionsTargetedLowPurity(double aPurity, double expectedPurityThreshold,
            @NotNull FindingStatus.Status expectedStatus,
            Set<FindingStatus.Issue> expectedErrors, Set<FindingStatus.Issue> expectedWarnings, int expectedResults)
    {
        FindingRecord original = findingRecordBuilder(aPurity, SequencingScope.WHOLE_GENOME)
                .metaProperties(TestFindingRecordFactory.metaPropertiesBuilder()
                        .sequencingScope(SequencingScope.TARGETED)
                        .build())
                .somaticGainDeletions(TestFindingFactory.driverFindingsBuilder(List.of(TestFindingFactory.gainDeletionBuilder()
                                        .somaticType(GainDeletion.Type.HOM_DEL)
                                        .build(),
                                TestFindingFactory.gainDeletionBuilder()
                                        .somaticType(GainDeletion.Type.GAIN)
                                        .build()))
                        .build())
                .build();
        FindingRecord converted = LowPurityTransformer.transform(original);

        DriverFindingList<?> gainDeletions = converted.somaticGainDeletions();
        assertFindingStatus(gainDeletions.status(), expectedStatus, expectedErrors, expectedWarnings);
        assertPurityThreshold(gainDeletions, expectedPurityThreshold);
        assertEquals(expectedResults, gainDeletions.findings().size());
    }

    private void assertQc(Qc qc)
    {
        assertFalse(qc.isPass());
        assertEquals(qc.errors(), new TreeSet<>(Set.of(Qc.QcStatus.LOW_PURITY)));
    }

    private void assertFindingItemOk(FindingItem<?> item)
    {
        assertFindingStatus(item.status(), FindingStatus.Status.OK, Set.of(), Set.of());
        assertPurityThreshold(item, null);
        assertNotNull(item.finding());
    }

    private void assertFindingItemLowPurity(FindingItem<?> item, @Nullable Double expectedPurityThreshold)
    {
        assertFindingItem(item, expectedPurityThreshold, Set.of(LOW_PURITY), Set.of());
    }

    private void assertFindingItem(FindingItem<?> item, @Nullable Double expectedPurityThreshold,
            Set<FindingStatus.Issue> expectedErrors, Set<FindingStatus.Issue> expectedWarnings)
    {
        assertFindingStatus(item.status(), FindingStatus.Status.NOT_RELIABLE, expectedErrors, expectedWarnings);
        assertPurityThreshold(item, expectedPurityThreshold);
        assertNull(item.finding());
    }

    private void assertFindingListOk(IFindingList<?> list)
    {
        assertFindingListOk(list, null, Set.of());
    }

    private void assertFindingListOk(IFindingList<?> list, @Nullable Double expectedPurityThreshold, Set<FindingStatus.Issue> warnings)
    {
        assertFindingStatus(list.status(), FindingStatus.Status.OK, Set.of(), warnings);
        assertPurityThreshold(list, expectedPurityThreshold);
        assertFalse(list.findings().isEmpty());
    }

    private void assertFindingListLowPurity(IFindingList<?> list, @Nullable Double expectedPurityThreshold)
    {
        assertFindingList(list, expectedPurityThreshold, Set.of(LOW_PURITY), Set.of());
    }

    private void assertFindingList(IFindingList<?> list, @Nullable Double expectedPurityThreshold,
            Set<FindingStatus.Issue> expectedErrors, Set<FindingStatus.Issue> expectedWarnings)
    {
        assertFindingStatus(list.status(), FindingStatus.Status.NOT_RELIABLE, expectedErrors, expectedWarnings);
        assertPurityThreshold(list, expectedPurityThreshold);
        assertTrue(list.findings().isEmpty());
    }

    private void assertPurityThreshold(@NotNull FindingItem<?> item, @Nullable Double expectedPurityThreshold)
    {
        assertPurityThreshold(expectedPurityThreshold, item.purityThreshold());
    }

    private void assertPurityThreshold(@NotNull IFindingList<?> list, @Nullable Double expectedPurityThreshold)
    {
        assertPurityThreshold(expectedPurityThreshold, list.purityThreshold());
    }

    private void assertPurityThreshold(@Nullable Double expectedPurityThreshold, @Nullable Double purityThreshold)
    {
        if(expectedPurityThreshold != null)
        {
            assertNotNull(purityThreshold);
            assertEquals(expectedPurityThreshold, purityThreshold, EPSILON);
        }
        else
        {
            assertNull(purityThreshold);
        }
    }

    private void assertFindingStatusOk(FindingStatus findingStatus)
    {
        assertFindingStatus(findingStatus, FindingStatus.Status.OK, Set.of(), Set.of());
    }

    private void assertFindingStatusLowPurity(FindingStatus findingStatus)
    {
        assertFindingStatus(findingStatus, FindingStatus.Status.NOT_RELIABLE, Set.of(LOW_PURITY), Set.of());
    }

    private void assertFindingStatus(FindingStatus findingStatus, @NotNull FindingStatus.Status expectedStatus,
            Set<FindingStatus.Issue> expectedErrors, Set<FindingStatus.Issue> expectedWarnings)
    {
        assertEquals(expectedStatus, findingStatus.status());
        assertEquals(expectedErrors, findingStatus.errors());
        assertEquals(expectedWarnings, findingStatus.warnings());
    }

    private static FindingRecordBuilder findingRecordBuilder(double purity, SequencingScope sequencingScope)
    {
        return TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .metaProperties(TestFindingRecordFactory.metaPropertiesBuilder()
                        .sequencingScope(sequencingScope)
                        .build())
                .qc(TestFindingFactory.qcBuilder()
                        .build())
                .purityPloidyFit(TestFindingFactory.purityPloidyFitBuilder()
                        .purity(purity)
                        .build());
    }
}
