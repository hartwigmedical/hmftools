package com.hartwig.hmftools.sage.common;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.IndexedBases;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefSequenceTest
{

    private static final ReferenceSequence REF_SEQ = new ReferenceSequence("1", 999, new byte[] { 'G', 'A', 'T', 'A', 'C', 'A' });

    @Test
    public void testRNA()
    {
        RefSequence rnaRefSequence = new RefSequence(REF_SEQ);
        SAMRecord samRecord = buildSamRecord(990, "20M", "GGGGGAAAAATTTTTCCCCC", "XXXXXXXXXXXXXXXXXXXX");
        IndexedBases indexedBases = rnaRefSequence.alignment();
        assertEquals(1000, indexedBases.Position);
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }
}
