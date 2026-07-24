package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.XA_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

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

        Read.AltAlignment firstAlt = read.altLoci().get(0);
        assertEquals("2", firstAlt.Region.Chromosome);
        assertEquals(5000, firstAlt.Region.start());
        assertEquals(5019, firstAlt.Region.end()); // 20M spans 20 reference bases
        assertFalse(firstAlt.Spliced);

        Read.AltAlignment secondAlt = read.altLoci().get(1);
        assertEquals("3", secondAlt.Region.Chromosome);
        assertEquals(8000, secondAlt.Region.start());
    }

    @Test
    public void testSplicedAltSpanFromCigar()
    {
        // an alt whose CIGAR has an N gap is flagged spliced and its reference span includes the skipped intron
        Read read = Read.from(createRead("5,+7000,10M100N10M,0;"));

        assertEquals(2, read.numLoci());

        Read.AltAlignment alt = read.altLoci().get(0);
        assertTrue(alt.Spliced);
        assertEquals(7000, alt.Region.start());
        assertEquals(7119, alt.Region.end()); // 10 + 100 + 10 = 120 reference bases
    }

    @Test
    public void testMalformedXaEntriesTolerated()
    {
        // a malformed and an empty entry are skipped, leaving the two valid alternate loci
        Read read = Read.from(createRead("2,+5000,20M,1;garbage;;4,+9000,20M,0;"));

        assertEquals(3, read.numLoci());
        assertEquals(2, read.altLoci().size());
        assertEquals("4", read.altLoci().get(1).Region.Chromosome);
        assertEquals(9000, read.altLoci().get(1).Region.start());
    }
}
