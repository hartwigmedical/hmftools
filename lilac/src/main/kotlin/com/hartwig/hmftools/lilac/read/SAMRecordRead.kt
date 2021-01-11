package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.lilac.ext.containsIndel
import com.hartwig.hmftools.lilac.sam.SamSlicer
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import java.util.function.Consumer
import kotlin.math.max
import kotlin.math.min


data class SAMRecordRead(val hlaNucIndexStart: Int, val readIndexStart: Int, val length: Int, val reverseStrand: Boolean, val samRecord: SAMRecord) : Read {

    private val hlaNucIndexEnd = hlaNucIndexStart + length - 1
    private val nucleotideIndices = IntRange(hlaNucIndexStart, hlaNucIndexEnd).toList()
    private val aminoAcidIndices =  AminoAcidIndices.indices(hlaNucIndexStart, hlaNucIndexEnd).toList()

    fun containsNucleotide(index: Int): Boolean {
        return index in nucleotideIndices
    }

    override fun nucleotideIndices(): Collection<Int> {
        return nucleotideIndices
    }

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidIndices.contains(index)
    }

    override fun aminoAcidIndices(): Collection<Int> {
        return aminoAcidIndices
    }

    override fun aminoAcid(index: Int, minQual: Int): Char {
        val startIndex = index * 3

        val firstBase = nucleotide(startIndex + 0, minQual)
        if (firstBase == '.') {
            return '.'
        }

        val secondBase = nucleotide(startIndex + 1, minQual)
        if (secondBase == '.') {
            return '.'
        }

        val thirdBase = nucleotide(startIndex + 2, minQual)
        if (thirdBase == '.') {
            return '.'
        }

        return Codons.aminoAcid(firstBase.toString() + secondBase + thirdBase)
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

    companion object {

        fun realign(hlaCodingRegionOffset: Int, region: GenomeRegion, reverseStrand: Boolean, bamFileName: String): List<SAMRecordRead> {
            val slicer = SamSlicer(1)
            val result = mutableListOf<SAMRecordRead>()
            SamReaderFactory.makeDefault().open(File(bamFileName)).use { samReader ->

                val consumer = Consumer<SAMRecord> { samRecord ->
                    if (!samRecord.containsIndel()) {
                        if (reverseStrand) {
                            result.add(realignReverseStrand(hlaCodingRegionOffset, region, samRecord))
                        } else {
                            result.add(realignForwardStrand(hlaCodingRegionOffset, region, samRecord))
                        }
                    }
                }


                slicer.slice(region, samReader, consumer)
            }
            return result
        }


        private fun realignForwardStrand(hlaExonOffset: Int, region: GenomeRegion, samRecord: SAMRecord): SAMRecordRead {
            val hlaExonStartPosition = region.start().toInt()
            val hlaExonEndPosition = region.end().toInt()

            val alignmentStart = samRecord.alignmentStart
            val alignmentEnd = samRecord.alignmentEnd

            val hlaStart = max(alignmentStart, hlaExonStartPosition)
            val hlaEnd = min(alignmentEnd, hlaExonEndPosition)
            val length = hlaEnd - hlaStart + 1

            val readIndex = samRecord.getReadPositionAtReferencePosition(hlaStart) - 1
            val hlaStartIndex = hlaStart - hlaExonStartPosition + hlaExonOffset

            return SAMRecordRead(hlaStartIndex, readIndex, length, false, samRecord)
        }

        private fun realignReverseStrand(hlaExonOffset: Int, region: GenomeRegion, samRecord: SAMRecord): SAMRecordRead {
            val hlaExonStartPosition = region.end().toInt()
            val hlaExonEndPosition = region.start().toInt()

            val alignmentStart = samRecord.alignmentStart
            val alignmentEnd = samRecord.alignmentEnd

            val hlaStart = min(alignmentEnd, hlaExonStartPosition)
            val hlaEnd = max(alignmentStart, hlaExonEndPosition)
            val length = hlaStart - hlaEnd + 1

            val readIndex = samRecord.getReadPositionAtReferencePosition(hlaStart) - 1
            val hlaStartIndex = hlaExonStartPosition - hlaStart + hlaExonOffset

            return SAMRecordRead(hlaStartIndex, readIndex, length, true, samRecord)
        }

    }

}
