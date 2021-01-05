package com.hartwig.hmftools.lilac.align

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.ext.containsIndel
import com.hartwig.hmftools.lilac.sam.SamSlicer
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import java.io.File
import java.util.function.Consumer
import kotlin.math.max
import kotlin.math.min


data class HlaAlign(val hlaIndexStart: Int, val readIndexStart: Int, val length: Int, val samRecord: SAMRecord) {
    val hlaIndexEnd = hlaIndexStart + length - 1


    fun containsIndex(index: Int): Boolean {
        return index >= hlaIndexStart && index <= hlaIndexStart + length
    }

    fun charAt(index: Int): Char {
        val adjustedIndex = index - hlaIndexStart + readIndexStart
        return samRecord.readBases[adjustedIndex].toChar()
    }

    companion object {

        fun realign(hlaCodingRegionOffset: Int, hlaCodingRegion: GenomeRegion, bamFileName: String): List<HlaAlign> {
            val slicer = SamSlicer(1)
            val result = mutableListOf<HlaAlign>()

            SamReaderFactory.makeDefault().open(File(bamFileName)).use { samReader ->

                val consumer = Consumer<SAMRecord> { samRecord ->
                    if (!samRecord.containsIndel()) {
                        result.add(realign(hlaCodingRegionOffset, hlaCodingRegion, samRecord))
                    }
                }

                slicer.slice(samReader, hlaCodingRegion, consumer)
            }

            return result
        }


        private fun realign(hlaCodingRegionOffset: Int, hlaCodingRegion: GenomeRegion, samRecord: SAMRecord): HlaAlign {
            val alignmentStart = samRecord.alignmentStart
            val alignmentEnd = samRecord.alignmentEnd

            val hlaStart = max(alignmentStart, hlaCodingRegion.start().toInt())
            val hlaEnd = min(alignmentEnd, hlaCodingRegion.end().toInt())
            val length = hlaEnd - hlaStart + 1

            val alignmentStartIndex = samRecord.getReadPositionAtReferencePosition(alignmentStart) - 1
            val hlaStartIndex = hlaStart - alignmentStart + alignmentStartIndex
            val hlaIndex = hlaStart - hlaCodingRegion.start().toInt() + hlaCodingRegionOffset

            return HlaAlign(hlaIndex, hlaStartIndex, length, samRecord)
        }

    }

}
