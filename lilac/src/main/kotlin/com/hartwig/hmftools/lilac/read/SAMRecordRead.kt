package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.bed.NamedBed
import com.hartwig.hmftools.common.genome.region.CodingRegions
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.lilac.LilacApplication2
import com.hartwig.hmftools.lilac.ext.containsIndel
import com.hartwig.hmftools.lilac.sam.SamSlicer
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import java.util.function.Consumer
import kotlin.math.max
import kotlin.math.min


data class SAMRecordRead(val gene: String, val hlaNucIndexStart: Int, val readIndexStart: Int, val length: Int, val reverseStrand: Boolean, val samRecord: SAMRecord) : Read {

    private val hlaNucIndexEnd = hlaNucIndexStart + length - 1
    private val nucleotideIndices = IntRange(hlaNucIndexStart, hlaNucIndexEnd).toList()
    private val aminoAcidIndices = AminoAcidIndices.indices(hlaNucIndexStart, hlaNucIndexEnd).toList()

    fun containsNucleotide(index: Int): Boolean {
        return index in nucleotideIndices
    }

    override fun nucleotideIndices(minQual: Int): Collection<Int> {
        return nucleotideIndices.filter { nucleotideQuality(it) >= minQual }
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

    fun nucleotideQuality(loci: Int): Int {
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

        fun readFromBam(transcript: HmfTranscriptRegion, bamFile: String): List<SAMRecordRead> {
            LilacApplication2.logger.info("... querying ${transcript.gene()} (${transcript.chromosome()}:${transcript.codingStart()}-${transcript.codingEnd()})")

            val reverseStrand = transcript.strand() == Strand.REVERSE
            val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)

            val realignedRegions = mutableListOf<SAMRecordRead>()
            var length = 0
            for (codingRegion in codingRegions) {
                realignedRegions.addAll(realign(transcript.gene(), length, codingRegion, reverseStrand, bamFile))
                length += codingRegion.bases().toInt()
//            println((length - 1) / 3)
            }
            return realignedRegions
        }

        private fun codingRegions(transcript: HmfTranscriptRegion): List<NamedBed> {
            return CodingRegions.codingRegions(transcript)
        }

        private fun realign(gene: String, hlaCodingRegionOffset: Int, region: GenomeRegion, reverseStrand: Boolean, bamFileName: String): List<SAMRecordRead> {
            val slicer = SamSlicer(1)
            val result = mutableListOf<SAMRecordRead>()
            SamReaderFactory.makeDefault().open(File(bamFileName)).use { samReader ->

                val consumer = Consumer<SAMRecord> { samRecord ->
                    if (!samRecord.containsIndel()) {
                        if (reverseStrand) {
                            result.add(realignReverseStrand(gene, hlaCodingRegionOffset, region, samRecord))
                        } else {
                            result.add(realignForwardStrand(gene, hlaCodingRegionOffset, region, samRecord))
                        }
                    }
                }


                slicer.slice(region, samReader, consumer)
            }
            return result
        }


        private fun realignForwardStrand(gene: String, hlaExonOffset: Int, region: GenomeRegion, samRecord: SAMRecord): SAMRecordRead {
            val hlaExonStartPosition = region.start().toInt()
            val hlaExonEndPosition = region.end().toInt()

            val alignmentStart = samRecord.alignmentStart
            val alignmentEnd = samRecord.alignmentEnd

            val hlaStart = max(alignmentStart, hlaExonStartPosition)
            val hlaEnd = min(alignmentEnd, hlaExonEndPosition)
            val length = hlaEnd - hlaStart + 1

            val readIndex = samRecord.getReadPositionAtReferencePosition(hlaStart) - 1
            val hlaStartIndex = hlaStart - hlaExonStartPosition + hlaExonOffset

            return SAMRecordRead(gene, hlaStartIndex, readIndex, length, false, samRecord)
        }

        private fun realignReverseStrand(gene: String, hlaExonOffset: Int, region: GenomeRegion, samRecord: SAMRecord): SAMRecordRead {
            val hlaExonStartPosition = region.end().toInt()
            val hlaExonEndPosition = region.start().toInt()

            val alignmentStart = samRecord.alignmentStart
            val alignmentEnd = samRecord.alignmentEnd

            val hlaStart = min(alignmentEnd, hlaExonStartPosition)
            val hlaEnd = max(alignmentStart, hlaExonEndPosition)
            val length = hlaStart - hlaEnd + 1

            val readIndex = samRecord.getReadPositionAtReferencePosition(hlaStart) - 1
            val hlaStartIndex = hlaExonStartPosition - hlaStart + hlaExonOffset

            return SAMRecordRead(gene, hlaStartIndex, readIndex, length, true, samRecord)
        }

    }

}
