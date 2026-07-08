package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.common.Read;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MultiMapReadTest
{
    private static final String READ_BASES = "ACGTACGTACGTACGTACGT"; // 20 bases, matches 20M

    private static SAMRecord createRead(final String xaTag)
    {
        SAMRecord record = createSamRecord(
                "READ_01", "1", 1000, READ_BASES, "20M", "1", 1200, false, false, null);

        if(xaTag != null)
            record.setAttribute(XA_ATTRIBUTE, xaTag);

        return record;
    }

    @Test
    public void testUniqueReadHasSingleLocus()
    {
        Read read = Read.from(createRead(null));
        assertEquals(1, read.numLoci());
        assertNull(read.altLoci());
    }

    @Test
    public void testXaParsedIntoAltLoci()
    {
        // standard lifted bwa XA: chr,+/-pos,CIGAR,NM; - the position sign encodes strand and must be stripped
        Read read = Read.from(createRead("2,+5000,20M,1;3,-8000,20M,2;"));

        assertEquals(3, read.numLoci());
        assertEquals(2, read.altLoci().size());

        ChrBaseRegion firstAlt = read.altLoci().get(0);
        assertEquals("2", firstAlt.Chromosome);
        assertEquals(5000, firstAlt.start());

        ChrBaseRegion secondAlt = read.altLoci().get(1);
        assertEquals("3", secondAlt.Chromosome);
        assertEquals(8000, secondAlt.start());
    }

    @Test
    public void testMalformedXaEntriesTolerated()
    {
        // a malformed and an empty entry are skipped, leaving the two valid alternate loci
        Read read = Read.from(createRead("2,+5000,20M,1;garbage;;4,+9000,20M,0;"));

        assertEquals(3, read.numLoci());
        assertEquals(2, read.altLoci().size());
        assertEquals("4", read.altLoci().get(1).Chromosome);
        assertEquals(9000, read.altLoci().get(1).start());
    }
}
