package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.isofox.FragmentAllocator.starFragmentCount;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.isofox.common.Read;

import org.junit.After;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class AlignerReadTest
{
    private static final double EPSILON = 1e-9;
    private static final String READ_BASES = "ACGTACGTACGTACGTACGT";

    private static SAMRecord createRead()
    {
        return createSamRecord("READ_01", "1", 1000, READ_BASES, "20M", "1", 1200, false, false, null);
    }

    @After
    public void resetAligner()
    {
        IsofoxConstants.STAR_ALIGNER = false;
    }

    @Test
    public void testStarFragmentWeightByMapQual()
    {
        assertEquals(1.0, starFragmentCount(255), EPSILON);
        assertEquals(1.0, starFragmentCount(4), EPSILON);
        assertEquals(0.5, starFragmentCount(3), EPSILON);
        assertEquals(0.33, starFragmentCount(2), EPSILON);
        assertEquals(0.2, starFragmentCount(1), EPSILON);
        assertEquals(0.1, starFragmentCount(0), EPSILON);
    }

    @Test
    public void testStarMultiMappedFromMapQual()
    {
        IsofoxConstants.STAR_ALIGNER = true;

        SAMRecord multi = createRead();
        multi.setMappingQuality(3);
        assertTrue(Read.from(multi).isMultiMapped());

        SAMRecord unique = createRead();
        unique.setMappingQuality(255);
        assertFalse(Read.from(unique).isMultiMapped());
    }

    @Test
    public void testBwaMultiMappedFromXa()
    {
        SAMRecord multi = createRead();
        multi.setAttribute(XA_ATTRIBUTE, "2,+5000,20M,1;");
        assertTrue(Read.from(multi).isMultiMapped());

        // a low-MAPQ read with no XA is not multi-mapped under bwa-tars (only XA drives it)
        SAMRecord unique = createRead();
        unique.setMappingQuality(0);
        assertFalse(Read.from(unique).isMultiMapped());
    }
}
