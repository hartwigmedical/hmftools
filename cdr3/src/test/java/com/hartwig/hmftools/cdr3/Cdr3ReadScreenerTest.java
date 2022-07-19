package com.hartwig.hmftools.cdr3;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class Cdr3ReadScreenerTest
{
    @Test
    public void testExtrapolateReadOffsetAtRefPosition()
    {
        // 50 bases alignment
        int alignmentStart = 1001;
        int alignmentEnd = 1050;

        // create a SAM record with soft clip on both sides
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        // soft clips on both sides
        record.setCigarString("30S50M20S");
        // 100 bases here
        record.setReadString("GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT");
        record.setReadNegativeStrandFlag(false);
        //record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);

        // make sure we set it up correctly
        assertEquals(31, record.getReadPositionAtReferencePosition(alignmentStart));
        assertEquals(80, record.getReadPositionAtReferencePosition(alignmentEnd));

        GenomeRegion mapped = record.getReadUnmappedFlag() ? null : GenomeRegions.create(record.getReferenceName(),
                record.getAlignmentStart(), record.getAlignmentEnd());

        int readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart, mapped);
        // start offset is 30, due to 30 left soft clip
        assertEquals(30, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd, mapped);
        assertEquals(79, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 30, mapped);
        // start offset is 0, due to 30 left soft clip
        assertEquals(0, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 20, mapped);
        // start offset is 99, due to 20 right soft clip
        assertEquals(99, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30, mapped);
        // out of range
        assertEquals(-1, readOffset);

        // test out of range
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 31, mapped);
        assertEquals(-1, readOffset);

        // this is over the end, so should get -1
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30, mapped);
        assertEquals(-1, readOffset);
    }
}
