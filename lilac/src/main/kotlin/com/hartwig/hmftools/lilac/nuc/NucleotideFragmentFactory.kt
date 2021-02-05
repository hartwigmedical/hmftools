package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.common.genome.bed.NamedBed
import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.sam.SAMCodingRecord
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import htsjdk.samtools.SAMRecord


class NucleotideFragmentFactory(inserts: List<HlaSequenceLoci>, deletes: List<HlaSequenceLoci>) {
    val insertSuffixTrees = inserts.map { Pair(it, SuffixTree(it.sequence())) }.toMap()
    val deleteSuffixTrees = deletes.map { Pair(it, SuffixTree(it.sequence())) }.toMap()


//    private val reverseStrand: Boolean, private val startLoci: Int, private val codingRegion: GenomeRegion

    fun doStuff(record: SAMRecord, reverseStrand: Boolean, codingRegionLoci: Int, codingRegion: NamedBed): NucleotideFragment? {
        val samCoding = SAMCodingRecord.create(codingRegion, record)
        val samCodingLength = samCoding.positionEnd - samCoding.positionStart + 1
        val samCodingStartLoci = if (reverseStrand) {
            codingRegionLoci + samCoding.positionEnd - codingRegion.end().toInt()
        } else {
            codingRegionLoci + samCoding.positionStart - codingRegion.start().toInt()
        }
        val samCodingEndLoci =  samCodingStartLoci + samCodingLength - 1



        if (record.readName == "A00260:132:H25FTDSXY:4:2143:22887:32894") {
            // INSERT
            println("sdf")
        }

        if (record.readName == "A00260:132:H25FTDSXY:4:2618:14190:27289") {
            // SOFT CLIP
            println("sdf")
        }

        if (samCoding.containsIndel() || samCoding.containsSoftClip()) {
            val allowedCodingStartLociRange = (samCodingStartLoci - samCoding.softClippedStart)..samCodingStartLoci

            var sequence = samCoding.read()
            if (reverseStrand) {
                sequence = sequence.map { it.reverseCompliment() }.joinToString ("").reversed()
            }

           val jon = insertSuffixTrees
                   .map { Pair(it.key, it.value.indices(sequence)) }
                   .filter { it.second.isNotEmpty() }
                   .filter { it.second.any {i -> allowedCodingStartLociRange.contains(i)}  }
            if (jon.isNotEmpty()) {
                println("sdf")
            }

        }



        return null
    }

    private fun Char.reverseCompliment(): Char {
        when (this) {
            'G' -> return 'C'
            'A' -> return 'T'
            'T' -> return 'A'
            'C' -> return 'G'
        }

        return this
    }


}