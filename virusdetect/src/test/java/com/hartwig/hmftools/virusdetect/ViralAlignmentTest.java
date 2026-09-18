package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;

public class ViralAlignmentTest
{
    private static final int CONTIG_LENGTH = 200;
    private static final ViralReference REFERENCE = reference();

    // from() pulls the tags, clips and aligned blocks out of the SAMRecord, and counts clipped bases toward divergence:
    // a 40-base soft clip plus one mismatch is a divergence of 41.
    @Test
    public void testExtractsFieldsAndCountsClipsInDivergence()
    {
        ViralAlignment alignment = ViralAlignment.from(record(30, "40S60M", 60, 1), REFERENCE);

        assertEquals("v1", alignment.contig().name());
        assertEquals(30, alignment.alignmentStart());
        assertEquals(89, alignment.alignmentEnd());
        assertEquals(40, alignment.leftClip());
        assertEquals(0, alignment.rightClip());
        assertEquals(60, alignment.alignerScore());
        assertEquals(41, alignment.divergence());   // NM 1 + 40 clipped bases
        assertEquals(1, alignment.alignedIntervals().size());
        assertEquals(30, alignment.alignedIntervals().get(0).referenceStart());
        assertEquals(60, alignment.alignedIntervals().get(0).length());
    }

    // A clip projecting past a contig boundary (beyond tolerance) straddles the circular origin; one staying inside does not.
    @Test
    public void testDetectsClipsOverContigEnds()
    {
        assertTrue(ViralAlignment.from(record(30, "40S60M", 60, 0), REFERENCE).clipsOverContigEnd());    // left clip projects to -10
        assertTrue(ViralAlignment.from(record(141, "60M40S", 60, 0), REFERENCE).clipsOverContigEnd());   // right clip projects to 240
        assertFalse(ViralAlignment.from(record(100, "40S60M", 60, 0), REFERENCE).clipsOverContigEnd());  // left clip projects to 60
    }

    @Test
    public void testThrowsWhenAlignmentScoreMissing()
    {
        SAMRecord record = record(30, "60M", null, 0);
        assertThrows(IllegalStateException.class, () -> ViralAlignment.from(record, REFERENCE));
    }

    @Test
    public void testThrowsWhenEditDistanceMissing()
    {
        SAMRecord record = record(30, "60M", 60, null);
        assertThrows(IllegalStateException.class, () -> ViralAlignment.from(record, REFERENCE));
    }

    private static SAMRecord record(int start, String cigar, @Nullable Integer alignerScore, @Nullable Integer editDistance)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("v1", CONTIG_LENGTH));

        SAMRecord record = new SAMRecord(header);
        record.setReadName("r");
        record.setReferenceName("v1");
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases("A".repeat(TextCigarCodec.decode(cigar).getReadLength()).getBytes());
        if(alignerScore != null)
        {
            record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, alignerScore.intValue());
        }
        if(editDistance != null)
        {
            record.setAttribute(NUM_MUTATONS_ATTRIBUTE, editDistance.intValue());
        }
        return record;
    }

    private static ViralReference reference()
    {
        List<ViralContig> contigs = List.of(new ViralContig("v1", CONTIG_LENGTH, "Virus 1", "Group 1"));
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(new SAMSequenceRecord("v1", CONTIG_LENGTH)));
        return new ViralReference(contigs, dictionary);
    }
}
