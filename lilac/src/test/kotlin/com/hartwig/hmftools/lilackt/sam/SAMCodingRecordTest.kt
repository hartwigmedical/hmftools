package com.hartwig.hmftools.lilackt.sam

import com.hartwig.hmftools.common.genome.region.GenomeRegions
import htsjdk.samtools.SAMRecord
import org.junit.Assert.assertEquals
import org.junit.Assert.assertTrue
import org.junit.Test

class SAMCodingRecordTest {

    private val shortCodingRegion = GenomeRegions.create("1", 1000, 1020)
    private val longCodingRegion = GenomeRegions.create("1", 1000, 2000)

    @Test
    fun testCodingRegionContainsMatchedRecord() {
        val samRecord = buildSamRecord(1100, "150M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1100, 1249, 0, 149)
    }

    @Test
    fun testLeftSoftClipWithin() {
        val samRecord = buildSamRecord(1100, "10S140M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 10, 0, 0, 0, 1090, 1239, 0, 149)
    }

    @Test
    fun testRightSoftClipWithin() {
        val samRecord = buildSamRecord(1100, "140M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 10, 0, 0, 1100, 1249, 0, 149)
    }

    @Test
    fun testRecordIsLeft() {
        val samRecord = buildSamRecord(900, "10S140M", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1000, 1039, 110, 149)
    }

    @Test
    fun testRecordIsLeftWithSomeRightSoftClip() {
        val samRecord = buildSamRecord(900, "10S130M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 10, 0, 0, 1000, 1039, 110, 149)
    }

    @Test
    fun testRecordIsRight() {
        val samRecord = buildSamRecord(1900, "140M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1900, 2000, 0, 100)
    }

    @Test
    fun testRecordIsRightWithSomeLeftSoftClip() {
        val samRecord = buildSamRecord(1910, "10S130M10S", 150)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 10, 0, 0, 0, 1900, 2000, 0, 100)
    }

    @Test
    fun testInsert() {
        val samRecord = buildSamRecord(1100, "50M3I50M", 103)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 3, 1100, 1199, 0, 102)
    }

    @Test
    fun testInsertOutsideCodingRegion() {
        val samRecord = buildSamRecord(1990, "50M3I50M", 103)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1990, 2000, 0, 10)
    }

    @Test
    fun testDelete() {
        val samRecord = buildSamRecord(1100, "50M3D50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, -3, 0, 1100, 1202, 0, 99)
        assertTrue(victim.containsIndel())
    }

    @Test
    fun testOppositeIndelsDoNotCancelOut() {
        val samRecord = buildSamRecord(1100, "50M3D1M3I50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertEquals(2, victim.indels.size)
        assertTrue(victim.containsIndel())
    }

    @Test
    fun testMultipleIndels() {
        val samRecord = buildSamRecord(1100, "50M3D1M4I50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertEquals(2, victim.indels.size)
        assertTrue(victim.containsIndel())
    }

    @Test
    fun testDeleteOutsideCodingRange() {
        val samRecord = buildSamRecord(1990, "50M3D50M", 100)
        val victim = SAMCodingRecord.create(false, longCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1990, 2000, 0, 10)
    }

    @Test
    fun testShortRegion() {
        val samRecord = buildSamRecord(900, "150M", 150)
        val victim = SAMCodingRecord.create(false, shortCodingRegion, samRecord)
        assertIntersect(victim, 0, 0, 0, 0, 1000, 1020, 100, 120)
    }


    fun assertIntersect(victim: SAMCodingRecord, softClippedStart: Int, softClippedEnd: Int, dels: Int, ins: Int, positionStart: Int, positionEnd: Int, readIndexStart: Int, readIndexEnd: Int) {
        assertEquals(softClippedStart, victim.softClippedStart)
        assertEquals(softClippedEnd, victim.softClippedEnd)
        assertEquals(dels, victim.indels.filter { it.isDelete }.map { it.length }.sum())
        assertEquals(ins, victim.indels.filter { it.isInsert }.map { it.length }.sum())
        assertEquals(positionStart, victim.positionStart)
        assertEquals(positionEnd, victim.positionEnd)
        assertEquals(readIndexStart, victim.readStart)
        assertEquals(readIndexEnd, victim.readEnd)
    }

    fun buildSamRecord(alignmentStart: Int, cigar: String, length: Int): SAMRecord {
        return buildSamRecord(alignmentStart, cigar,
                (0 until length).map { 'G' }.joinToString(""),
                (0 until length).map { '#' }.joinToString(""))
    }

    fun buildSamRecord(alignmentStart: Int, cigar: String, readString: String, qualities: String): SAMRecord {
        val record = SAMRecord(null)
        record.readName = "READNAME"
        record.alignmentStart = alignmentStart
        record.cigarString = cigar
        record.readString = readString
        record.readNegativeStrandFlag = false
        record.baseQualityString = qualities
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.readPairedFlag = true
        return record
    }

}