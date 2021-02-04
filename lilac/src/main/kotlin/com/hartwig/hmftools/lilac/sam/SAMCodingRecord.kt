package com.hartwig.hmftools.lilac.sam

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import kotlin.math.max
import kotlin.math.min

data class SAMCodingRecord(
        val softClipped: Int, val deleted: Int, val inserted: Int,
        val positionStart: Int, val positionEnd: Int,
        val readStart: Int, val readEnd: Int,
        val record: SAMRecord) {

    companion object {
        fun create(codingRegion: GenomeRegion, record: SAMRecord): SAMCodingRecord {
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
            if (positionStart == alignmentStart && softClipStart > 0) {
                val earliestStart = max(codingRegion.start().toInt(), recordStart)
                readIndexStart = readIndexStart - positionStart + earliestStart
                positionStart = earliestStart
            }

            // Add soft clip end
            if (positionEnd == alignmentEnd && softClipEnd > 0) {
                val latestEnd = min(codingRegion.end().toInt(), recordEnd)
                readIndexEnd = readIndexEnd + latestEnd - positionEnd
                positionEnd = latestEnd
            }

            val softClipped = max(0, positionEnd - alignmentEnd) + max(alignmentStart - positionStart, 0)

            return SAMCodingRecord(softClipped, record.deletes(), record.inserts(), positionStart, positionEnd, readIndexStart, readIndexEnd, record)
        }

        private fun SAMRecord.softClipStart(): Int {
            return if (this.cigar.firstCigarElement.operator == CigarOperator.S) this.cigar.firstCigarElement.length else 0
        }

        private fun SAMRecord.softClipEnd(): Int {
            return if (this.cigar.lastCigarElement.operator == CigarOperator.S) this.cigar.lastCigarElement.length else 0
        }

        private fun SAMRecord.deletes(): Int {
            return this.cigar.cigarElements.filter { it.operator == CigarOperator.D }.map { it.length }.sum()
        }

        private fun SAMRecord.inserts(): Int {
            return this.cigar.cigarElements.filter { it.operator == CigarOperator.I }.map { it.length }.sum()
        }
    }

}