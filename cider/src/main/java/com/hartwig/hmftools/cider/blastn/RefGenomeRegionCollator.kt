package com.hartwig.hmftools.cider.blastn

import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.logging.log4j.LogManager
import java.lang.Math.max
import java.lang.Math.min
import java.util.*
import kotlin.Comparator

data class RefGenomeRegion(val contig: String, var start: Int, var endExclusive: Int, var sequence: String)
{
    init
    {
        require(endExclusive - start == sequence.length)
    }
}

class RefGenomeRegionCollator
{
    enum class Column
    {
        contig,
        start,
        endExclusive,
        sequence,
        count
    }

    data class RefGenomeRegionCount(val refGenomeRegion: RefGenomeRegion, var count: Int)

    val refGenomeRegions: MutableMap<String, MutableList<RefGenomeRegionCount>> = HashMap()

    fun addRefGenomeRegion(regionToAdd: RefGenomeRegion)
    {
        // see if regions already exist and also try to merge
        val regionList: MutableList<RefGenomeRegionCount> = refGenomeRegions.getOrPut(regionToAdd.contig, { ArrayList() })
        val regionCountToAdd = RefGenomeRegionCount(regionToAdd, 1)

        var i = Collections.binarySearch(regionList, regionCountToAdd, regionStartComparator)

        if (i < 0)
        {
            i = -(i + 1)
        }

        var insertNew = true

        // We want to check the regions before and after the insertion point and see if they overlap
        // this way we avoid unnecessary costly arraylist insert
        if (i > 0)
        {
            // we found the spot where the regionToAdd.start is larger than the regionList[i].start
            // but it is possible that the previous region overlaps with this region, we must check this
            val regionCount = regionList[i - 1]
            val region = regionCount.refGenomeRegion

            // condition from binary search
            require(regionToAdd.start > region.start)

            if (regionToAdd.start <= region.endExclusive)
            {
                mergeRegionCounts(regionCountToAdd, regionCount)
                i--
                insertNew = false
            }
        }
        if (insertNew && i < regionList.size)
        {
            // check if overlap
            val regionCount = regionList[i]
            val region = regionCount.refGenomeRegion

            // condition from binary search
            require(regionToAdd.start <= region.start)

            if (regionToAdd.endExclusive >= region.start)
            {
                // overlap
                mergeRegionCounts(regionCountToAdd, regionCount)
                insertNew = false
            }
        }

        if (insertNew)
        {
            regionList.add(i, regionCountToAdd)
        }

        val regionCount = regionList[i]
        val region = regionCount.refGenomeRegion

        // we might be able to swallow up later regions as well
        val j = i + 1
        while (j < regionList.size)
        {
            val regionCount1 = regionList[j]
            val region1 = regionCount1.refGenomeRegion
            if (region1.start <= region.endExclusive)
            {
                mergeRegionCounts(regionCount1, regionCount)
                regionList.removeAt(j)
            }
            else
            {
                break
            }
        }
    }

    fun readFromTsv(tsvPath: String)
    {
        FileWriterUtils.createBufferedReader(tsvPath).use { reader ->

            val csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                .build()

            for (record in csvFormat.parse(reader))
            {
                val region = RefGenomeRegion(
                    contig = record[Column.contig],
                    start = record[Column.start].toInt(),
                    endExclusive = record[Column.endExclusive].toInt(),
                    sequence = record[Column.sequence]
                )

                val count = record[Column.count].toInt()

                val l = refGenomeRegions.getOrPut(region.contig, { ArrayList<RefGenomeRegionCount>() })
                l.add(RefGenomeRegionCount(region, count))
            }
        }

        sLogger.info("loaded {} regions from {}",
            refGenomeRegions.values.sumOf { l -> l.size },
            tsvPath)
    }

    fun writeToTsv(tsvPath: String)
    {
        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        csvFormat.print(FileWriterUtils.createBufferedWriter(tsvPath)).use { printer: CSVPrinter ->
            for ((_, regionList) in refGenomeRegions.entries)
            {
                for (regionCount in regionList)
                {
                    val region = regionCount.refGenomeRegion
                    for (c in Column.values())
                    {
                        when (c)
                        {
                            Column.contig -> printer.print(region.contig)
                            Column.start -> printer.print(region.start)
                            Column.endExclusive -> printer.print(region.endExclusive)
                            Column.sequence -> printer.print(region.sequence)
                            Column.count -> printer.print(regionCount.count)
                        }
                    }
                    printer.println()
                }
            }
        }
    }

    companion object
    {
        val sLogger = LogManager.getLogger(RefGenomeRegionCollator::class.java)
        val regionStartComparator: Comparator<RefGenomeRegionCount> = Comparator.comparingInt {
                obj: RefGenomeRegionCount -> obj.refGenomeRegion.start }

        private fun mergeRegionCounts(
            regionCountFrom: RefGenomeRegionCount,
            regionCountTo: RefGenomeRegionCount
        )
        {
            // these regions overlap
            val regionFrom: RefGenomeRegion = regionCountFrom.refGenomeRegion
            val regionTo: RefGenomeRegion = regionCountTo.refGenomeRegion

            // check that the overlap regions are identical
            val overlapSeq1 = regionFrom.sequence.substring(
                max(regionTo.start - regionFrom.start, 0),
                regionFrom.sequence.length - max(regionFrom.endExclusive - regionTo.endExclusive, 0)
            )

            val overlapSeq2 = regionTo.sequence.substring(
                max(regionFrom.start - regionTo.start, 0),
                regionTo.sequence.length - max(regionTo.endExclusive - regionFrom.endExclusive, 0)
            )

            if (overlapSeq1 != overlapSeq2)
            {
                // test again without the non standard DNA codes
                val regex = Regex("[^ACTG]")
                if (overlapSeq1.replace(regex, "") != overlapSeq2.replace(regex, ""))
                {
                    sLogger.fatal("{} != {}", overlapSeq1, overlapSeq2)
                    throw RuntimeException()
                }
            }

            regionTo.sequence = regionFrom.sequence.take(max(regionTo.start - regionFrom.start, 0)) +
                    regionTo.sequence +
                    regionFrom.sequence.takeLast(max(regionFrom.endExclusive - regionTo.endExclusive, 0))
            regionTo.start = min(regionFrom.start, regionTo.start)
            regionTo.endExclusive = max(regionFrom.endExclusive, regionTo.endExclusive)

            regionCountTo.count += regionCountFrom.count
        }
    }
}