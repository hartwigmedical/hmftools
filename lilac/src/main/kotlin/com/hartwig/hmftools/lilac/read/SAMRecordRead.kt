package com.hartwig.hmftools.lilac.read

import htsjdk.samtools.SAMRecord


data class SAMRecordRead(val gene: String, val hlaNucIndexStart: Int, val readIndexStart: Int, val length: Int, val reverseStrand: Boolean, val samRecord: SAMRecord) : Read {

    private val hlaNucIndexEnd = hlaNucIndexStart + length - 1
    private val nucleotideIndices = IntRange(hlaNucIndexStart, hlaNucIndexEnd).toList()

    override fun nucleotideIndices(minQual: Int): Collection<Int> {
        return nucleotideIndices.filter { nucleotideQuality(it) >= minQual }
    }

    override fun nucleotideIndices(): Collection<Int> {
        return nucleotideIndices
    }

    private fun nucleotideQuality(loci: Int): Int {
        val readIndex = readIndex(loci)
        return samRecord.baseQualities[readIndex].toInt()
    }

    private fun readIndex(loci: Int): Int {
        return if (reverseStrand) readIndexStart - loci + hlaNucIndexStart else loci - hlaNucIndexStart + readIndexStart
    }

    override fun nucleotide(loci: Int): Char {
        val readIndex = readIndex(loci)
        val readBase = samRecord.readBases[readIndex].toChar()
        return if (reverseStrand) readBase.reverseCompliment() else readBase
    }

    fun quality(loci: Int): Int {
        val readIndex = readIndex(loci)
        return samRecord.baseQualities[readIndex].toInt()
    }

    override fun nucleotide(index: Int, minQual: Int): Char {
        val adjustedIndex = if (reverseStrand) readIndexStart - index + hlaNucIndexStart else index - hlaNucIndexStart + readIndexStart
        val quality = samRecord.baseQualities[adjustedIndex]
        if (quality < minQual) {
            return '.'
        }
        val readBase = samRecord.readBases[adjustedIndex].toChar()
        return if (reverseStrand) readBase.reverseCompliment() else readBase
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
