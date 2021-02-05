package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import htsjdk.samtools.SAMRecord
import org.junit.Test

class NucleotideFragmentFactoryTest {

    private val reference = HlaSequence("A*01:01", "GATTTACA..CAT..ATC")
    private val delete = HlaSequence("A*01:02", "--...---..---..---")
    private val insert = HlaSequence("A*01:03", "--------AA---TT---")


    @Test
    fun testStuff() {
        val hlaSequences = HlaSequenceLoci.create(listOf(reference, delete, insert))
        val hlaSequenceWithInserts = hlaSequences.filter { it.containsInserts() }
        val hlaSequenceWithDeletes = hlaSequences.filter { it.containsDeletes() }


    }

    fun createSAMRecord(position: Int) {

//        SAMRecord()

    }

    fun buildSamRecord(alignmentStart: Int, cigar: String, readString: String): SAMRecord {
        return buildSamRecord(alignmentStart, cigar, readString, readString.map { '#' }.joinToString(""))
    }

    fun buildSamRecord(alignmentStart: Int, cigar: String, readString: String, qualities: String): SAMRecord {
        val record = SAMRecord(null)
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