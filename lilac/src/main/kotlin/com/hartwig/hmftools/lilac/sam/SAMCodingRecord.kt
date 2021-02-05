package com.hartwig.hmftools.lilac.sam

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.samtools.CigarHandler
import com.hartwig.hmftools.common.samtools.CigarTraversal
import htsjdk.samtools.CigarElement
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord
import kotlin.math.max
import kotlin.math.min

class JonJon : CigarHandler {

}

data class SAMCodingRecord(
        val softClippedStart: Int, val softClippedEnd: Int,
        val deleted: Int, val inserted: Int,
        val positionStart: Int, val positionEnd: Int,
        val readStart: Int, val readEnd: Int,
        val record: SAMRecord) {

    fun containsSoftClip(): Boolean {
        return softClippedStart > 0 || softClippedEnd > 0
    }

    fun containsIndel(): Boolean {
        return deleted > 0 || inserted > 0
    }

    fun read(): String {
        val forwardBuilder = StringBuilder()
        for (i in readStart..readEnd) {
            forwardBuilder.append(record.readBases[i].toChar())
        }

        return forwardBuilder.toString()
    }

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

            val softClippedStart = max(alignmentStart - positionStart, 0)
            val softClippedEnd = max(0, positionEnd - alignmentEnd)
            val (insertCount, deleteCount) = indels(positionStart, positionEnd, record)

            return SAMCodingRecord(softClippedStart, softClippedEnd, deleteCount, insertCount, positionStart, positionEnd, readIndexStart, readIndexEnd, record)
        }

        private fun indels(startPosition: Int, endPosition: Int, record: SAMRecord): Pair<Int, Int> {
            var insertCount = 0
            var deleteCount = 0

            val handler = object : CigarHandler {

                override fun handleInsert(record: SAMRecord, element: CigarElement, readIndex: Int, refPosition: Int) {
                    if (refPosition in startPosition..endPosition) {
                        insertCount += element.length
                    }
                }

                override fun handleDelete(record: SAMRecord, element: CigarElement, readIndex: Int, refPosition: Int) {
                    if (refPosition in startPosition..endPosition) {
                        deleteCount += element.length
                    }
                }
            }

            CigarTraversal.traverseCigar(record, handler)
            return Pair(insertCount, deleteCount)
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