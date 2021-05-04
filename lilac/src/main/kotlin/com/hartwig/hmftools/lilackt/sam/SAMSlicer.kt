package com.hartwig.hmftools.lilackt.sam

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReader
import java.util.function.Consumer

class SAMSlicer(private val minMappingQuality: Int) {

    fun slice(region: GenomeRegion, samReader: SamReader, consumer: Consumer<SAMRecord>) {
        slice(region.chromosome(), region.start().toInt(), region.end().toInt(), samReader, consumer)
    }

    fun slice(contig: String, start: Int, end: Int, samReader: SamReader, consumer: Consumer<SAMRecord>) {
        val iterator = samReader.query(contig, start, end, false)
        for (samRecord in iterator) {
            if (samRecord.meetsQualityRequirements()) {
                consumer.accept(samRecord)
            }
        }
    }

    fun queryMates(samReader: SamReader, records: List<SAMRecord>): List<SAMRecord> {
        return records
                .mapNotNull { samReader.queryMate(it) }
                .filter { it.meetsQualityRequirements() }
    }

    private fun SAMRecord.meetsQualityRequirements(): Boolean {
        return this.mappingQuality >= minMappingQuality
                && !this.readUnmappedFlag
                && !this.duplicateReadFlag
                && !this.isSecondaryOrSupplementary
    }

}