package com.hartwig.hmftools.lilackt.sam

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.samtools.CigarHandler
import com.hartwig.hmftools.common.samtools.CigarTraversal
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min

data class SAMCodingRecord(
        val id: String,
        val softClippedStart: Int, val softClippedEnd: Int,
        val indels: List<Indel>,
        val positionStart: Int, val positionEnd: Int,
        val readStart: Int, val readEnd: Int,
        val record: SAMRecord,
        val reverseStrand: Boolean) {

    fun maxIndelSize(): Int {
        return indels.map { abs(it.length) }.max() ?: 0
    }

    fun containsSoftClip(): Boolean {
        return softClippedStart > 0 || softClippedEnd > 0
    }

    fun containsIndel(): Boolean {
        return indels.isNotEmpty()
    }

    fun codingRegionRead(reverseCompliment: Boolean): CharArray {
        return if (reverseCompliment)
            forwardRead().map { it.reverseCompliment() }.reversed().toCharArray()
        else
            forwardRead()
    }

    fun codingRegionQuality(reverseCompliment: Boolean): IntArray {
        return if (reverseCompliment)
            forwardQuality().reversed().toIntArray()
        else
            forwardQuality()
    }

    private fun forwardRead(): CharArray {
        return (readStart..readEnd).map { record.readBases[it].toChar() }.toCharArray()
    }

    private fun forwardQuality(): IntArray {
        return (readStart..readEnd).map { record.baseQualities[it].toInt() }.toIntArray()
    }

    fun alignmentsOnly(): List<SAMCodingRecord> {
        val outerRegion = GenomeRegions.create(record.contig, positionStart.toLong(), positionEnd.toLong())
        val result = record.alignmentBlocks
                .map { GenomeRegions.create(record.contig, it.referenceStart.toLong(), it.referenceStart.toLong() + it.length - 1) }
                .filter { outerRegion.overlaps(it) }
                .map { GenomeRegions.create(record.contig, max(it.start(), outerRegion.start()), min(it.end(), outerRegion.end())) }
                .map { create(this.reverseStrand, it, record, false, false) }

        return result

    }

    companion object {


        fun create(reverseStrand: Boolean, codingRegion: GenomeRegion, record: SAMRecord, includeSoftClips: Boolean = true, includeIndels: Boolean = true): SAMCodingRecord {

            val softClipStart = record.softClipStart()
            val softClipEnd = record.softClipEnd()

            val alignmentStart = record.alignmentStart
            val alignmentEnd = record.alignmentEnd

            val recordStart = alignmentStart - softClipStart
            val recordEnd = alignmentEnd + softClipEnd

            // Limit positions to coding regions
            var positionStart = max(codingRegion.start().toInt(), alignmentStart)
            var positionEnd = min(codingRegion.end().toInt(), alignmentEnd)
            var readIndexStart = record.getReadPositionAtReferencePosition(positionStart, true) - 1
            var readIndexEnd = record.getReadPositionAtReferencePosition(positionEnd, true) - 1
            positionStart = record.getReferencePositionAtReadPosition(readIndexStart + 1)
            positionEnd = record.getReferencePositionAtReadPosition(readIndexEnd + 1)

            // Add soft clip start
            if (positionStart == alignmentStart && softClipStart > 0 && includeSoftClips) {
                val earliestStart = max(codingRegion.start().toInt(), recordStart)
                readIndexStart = readIndexStart - positionStart + earliestStart
                positionStart = earliestStart
            }

            // Add soft clip end
            if (positionEnd == alignmentEnd && softClipEnd > 0 && includeSoftClips) {
                val latestEnd = min(codingRegion.end().toInt(), recordEnd)
                readIndexEnd = readIndexEnd + latestEnd - positionEnd
                positionEnd = latestEnd
            }

            val softClippedStart = max(alignmentStart - positionStart, 0)
            val softClippedEnd = max(0, positionEnd - alignmentEnd)
            val indels = if (includeIndels) indels(positionStart, positionEnd, record) else listOf()

            return SAMCodingRecord(record.readName, softClippedStart, softClippedEnd, indels, positionStart, positionEnd, readIndexStart, readIndexEnd, record, reverseStrand)
        }

        private fun indels(startPosition: Int, endPosition: Int, record: SAMRecord): List<Indel> {
            val indels = mutableListOf<Indel>()

            val handler = object : CigarHandler {

                override fun handleInsert(record: SAMRecord, element: CigarElement, readIndex: Int, refPosition: Int) {
                    if (refPosition in startPosition..endPosition) {
                        val base = record.readBases[readIndex].toChar()
                        val insert = Indel(record.contig, refPosition, base.toString(), record.readString.substring(readIndex, readIndex + element.length + 1))
                        indels.add(insert)
                    }
                }

                override fun handleDelete(record: SAMRecord, element: CigarElement, readIndex: Int, refPosition: Int) {
                    if (refPosition in startPosition..endPosition) {
                        val base = record.readBases[readIndex].toChar()
                        val delete = Indel(record.contig, refPosition, base + "N".repeat(element.length), base.toString())
                        indels.add(delete)
                    }
                }
            }

            CigarTraversal.traverseCigar(record, handler)
            return indels
        }

        private fun SAMRecord.softClipStart(): Int {
            return if (this.cigar.firstCigarElement.operator == CigarOperator.S) this.cigar.firstCigarElement.length else 0
        }

        private fun SAMRecord.softClipEnd(): Int {
            return if (this.cigar.lastCigarElement.operator == CigarOperator.S) this.cigar.lastCigarElement.length else 0
        }

        fun Char.reverseCompliment(): Char {
            when (this) {
                'G' -> return 'C'
                'A' -> return 'T'
                'T' -> return 'A'
                'C' -> return 'G'
            }

            return this
        }
    }

}