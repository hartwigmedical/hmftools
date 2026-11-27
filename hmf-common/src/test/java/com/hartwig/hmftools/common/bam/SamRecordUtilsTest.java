package com.hartwig.hmftools.common.bam;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.calculateNmTag;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.setNmTagIfMissing;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecordUnpaired;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SamRecordUtilsTest
{
    private static final String REF_BASES = "A".repeat(1_100);

    @Test
    public void testCalculateNmTagAcrossCigarOperations()
    {
        TrackingRefGenome refGenome = new TrackingRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES);

        String readBases = "AACAGGAAAAAAAACAACC";
        SAMRecord record = createSamRecordUnpaired(
                "read", CHR_1, 10, readBases, "4M2I3M1D4M1000N4M2S", false, false, null);

        assertEquals(5, calculateNmTag(record, refGenome));
        assertEquals(15, refGenome.RequestedBaseCount);

        assertTrue(setNmTagIfMissing(record, refGenome));
        assertEquals(5, record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());

        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, 42);
        assertFalse(setNmTagIfMissing(record, refGenome));
        assertEquals(42, record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
    }

    @Test
    public void testZeroBasedRefGenome()
    {
        MockRefGenome refGenome = new MockRefGenome(false);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES);

        SAMRecord record = createSamRecordUnpaired("read", CHR_1, 10, "AAAA", "4M", false, false, null);

        assertTrue(setNmTagIfMissing(record, refGenome));
        assertEquals(0, record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE).intValue());
    }

    @Test
    public void testInvalidRecordsAreIgnored()
    {
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES);

        SAMRecord unmappedRecord = createSamRecordUnpaired("unmapped", "*", 0, "AAAA", "*", false, false, null);
        assertFalse(setNmTagIfMissing(unmappedRecord, refGenome));
        assertNull(unmappedRecord.getAttribute(NUM_MUTATONS_ATTRIBUTE));

        SAMRecord noBasesRecord = createSamRecordUnpaired("noBases", CHR_1, 10, "", "*", false, false, null);
        assertFalse(setNmTagIfMissing(noBasesRecord, refGenome));
        assertNull(noBasesRecord.getAttribute(NUM_MUTATONS_ATTRIBUTE));
    }

    private static class TrackingRefGenome extends MockRefGenome
    {
        public int RequestedBaseCount;

        public TrackingRefGenome(boolean oneBasedIndexing)
        {
            super(oneBasedIndexing);
            RequestedBaseCount = 0;
        }

        @Override
        public byte[] getBases(final String chromosome, int posStart, int posEnd)
        {
            RequestedBaseCount += posEnd - posStart + 1;
            return super.getBases(chromosome, posStart, posEnd);
        }
    }
}
