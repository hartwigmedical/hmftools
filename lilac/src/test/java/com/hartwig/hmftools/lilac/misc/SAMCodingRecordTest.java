package com.hartwig.hmftools.lilac.misc;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.lilac.read.SAMCodingRecord;

import htsjdk.samtools.SAMRecord;

public class SAMCodingRecordTest
{
    private GenomeRegion shortCodingRegion = GenomeRegions.create("1", 1000, 1020);
    private GenomeRegion longCodingRegion = GenomeRegions.create("1", 1000, 2000);

            /*
    @Test
    public void testCodingRegionContainsMatchedRecord()
    {
        val samRecord = buildSamRecord(1100, "150M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1100, 1249, 0, 149)
    }

    @Test
    public void testLeftSoftClipWithin()
    {
        val samRecord = buildSamRecord(1100, "10S140M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 10, 0, 0, 0, 1090, 1239, 0, 149)
    }

    @Test
    public void testRightSoftClipWithin()
    {
        val samRecord = buildSamRecord(1100, "140M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 10, 0, 0, 1100, 1249, 0, 149)
    }

    @Test
    public void testRecordIsLeft()
    {
        val samRecord = buildSamRecord(900, "10S140M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1000, 1039, 110, 149)
    }

    @Test
    public void testRecordIsLeftWithSomeRightSoftClip()
    {
        val samRecord = buildSamRecord(900, "10S130M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 10, 0, 0, 1000, 1039, 110, 149)
    }

    @Test
    public void testRecordIsRight()
    {
        val samRecord = buildSamRecord(1900, "140M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1900, 2000, 0, 100)
    }

    @Test
    public void testRecordIsRightWithSomeLeftSoftClip()
    {
        val samRecord = buildSamRecord(1910, "10S130M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 10, 0, 0, 0, 1900, 2000, 0, 100)
    }

    @Test
    public void testInsert()
    {
        val samRecord = buildSamRecord(1100, "50M3I50M", 103)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 3, 1100, 1199, 0, 102)
    }

    @Test
    public void testInsertOutsideCodingRegion()
    {
        val samRecord = buildSamRecord(1990, "50M3I50M", 103)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1990, 2000, 0, 10)
    }

    @Test
    public void testDelete()
    {
        val samRecord = buildSamRecord(1100, "50M3D50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, -3, 0, 1100, 1202, 0, 99)
        assertTrue(victim.containsIndel())
    }

    @Test
    public void testOppositeIndelsDoNotCancelOut()
    {
        val samRecord = buildSamRecord(1100, "50M3D1M3I50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertEquals(2, victim.indels.size)
        assertTrue(victim.containsIndel())
    }

    @Test
    public void testMultipleIndels()
    {
        val samRecord = buildSamRecord(1100, "50M3D1M4I50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertEquals(2, victim.indels.size)
        assertTrue(victim.containsIndel())
    }

    @Test
    public void testDeleteOutsideCodingRange()
    {
        val samRecord = buildSamRecord(1990, "50M3D50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1990, 2000, 0, 10)
    }

    @Test
    public void testShortRegion()
    {
        val samRecord = buildSamRecord(900, "150M", 150)
        val victim = SAMCodingRecord.create(false, shortCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1000, 1020, 100, 120)
    }

    private void assertIntersect(
            SAMCodingRecord victim, int softClippedStart, int softClippedEnd, int dels, int ins,
            int positionStart, int positionEnd, int readIndexStart, int readIndexEnd)
    {
        assertEquals(softClippedStart, victim.SoftClippedStart);
        assertEquals(softClippedEnd, victim.SoftClippedEnd);

        assertEquals(dels, victim.getIndels().stream().filter(x -> x.IsDelete).mapToInt(x -> x.Length).sum();
        assertEquals(ins, victim.getIndels().stream().filter(x -> x.IsInsert).mapToInt(x -> x.Length).sum();

        assertEquals(positionStart, victim.PositionStart);
        assertEquals(positionEnd, victim.PositionEnd);
        assertEquals(readIndexStart, victim.ReadStart);
        assertEquals(readIndexEnd, victim.ReadEnd);
    }

             */

    private SAMRecord buildSamRecord(int alignmentStart, String cigar, int length)
    {
        //(0until length).map {'G'}.joinToString(""), (0 until length).map {'#'}.joinToString(""))

        String readStr = "";
        String qualsStr = "";

        for(int i = 0; i < length; ++i)
        {
            readStr += "G";
            qualsStr += "#";

        }

        return buildSamRecord(alignmentStart, cigar, readStr, qualsStr);
    }

    private SAMRecord buildSamRecord(int alignmentStart, String cigar, String readString, String qualities)
    {
        SAMRecord record = new SAMRecord(null);
        record.setReadName("READNAME");
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        return record;
    }

}
