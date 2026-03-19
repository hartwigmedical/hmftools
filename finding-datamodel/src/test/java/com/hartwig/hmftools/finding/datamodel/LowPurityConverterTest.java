package com.hartwig.hmftools.finding.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.finding.util.LowPurityConverter;

import org.junit.Test;

public class LowPurityConverterTest
{
    @Test
    public void testConvert()
    {
        FindingRecord original = TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .qc(TestFindingFactory.qcBuilder()
                        .status(Set.of(Qc.QCStatus.PASS, Qc.QCStatus.WARN_LOW_PURITY))
                        .build())
                .build();
        FindingRecord converted = LowPurityConverter.convert(original);

        assertFindingList(converted.somaticDisruptions());
        assertFindingList(converted.somaticGainDeletions());
        assertFindingList(converted.viruses());

        assertFindingItem(converted.homologousRecombination());
        assertFindingItem(converted.microsatelliteStability());
        assertFindingItem(converted.tumorMutationalBurden());
        assertFindingItem(converted.tumorMutationalLoad());
    }

    private void assertFindingItem(FindingItem<?> item)
    {
        assertEquals(FindingsStatus.NOT_RELIABLE, item.status());
        assertNull(item.finding());
    }

    private void assertFindingList(DriverFindingList<?> item)
    {
        assertEquals(FindingsStatus.NOT_RELIABLE, item.status());
        assertTrue(item.findings().isEmpty());
    }
}
