package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

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
