package com.hartwig.hmftools.cobalt.targeted

import com.hartwig.hmftools.common.genome.position.GenomePosition
import com.hartwig.hmftools.common.genome.position.GenomePositions
import com.hartwig.hmftools.common.utils.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import java.io.IOException
import java.util.*


class TargetRegionEnrichment
{
    private val mTargetedRegions: MutableList<GenomePosition> = ArrayList()
    private val mTargetRelativeEnrichment: MutableMap<GenomePosition, Double> = TreeMap(GenomePosition::compare)
    val targetedRegions: List<GenomePosition>
        get() = mTargetedRegions
    val targetRelativeEnrichment: Map<GenomePosition, Double>
        get() = mTargetRelativeEnrichment

    companion object
    {
        private const val DELIMITER = '\t'

        @JvmStatic
        @Throws(IOException::class)
        fun fromTsv(fileName: String): TargetRegionEnrichment
        {
            val targetRegionEnrichment = TargetRegionEnrichment()

            FileWriterUtils.createBufferedReader(fileName).use { reader ->

                val format: CSVFormat = CSVFormat.Builder.create()
                    .setDelimiter(DELIMITER)
                    .setRecordSeparator('\n')
                    .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                    .build()
                val records: Iterable<CSVRecord> = format.parse(reader)

                for (record: CSVRecord in records)
                {
                    val chromosome: String = record["chromosome"].intern()
                    val position: Int = record["position"].toDouble().toInt()
                    val relativeEnrichment: Double? = record["relativeEnrichment"].toDoubleOrNull()
                    val genomePosition = GenomePositions.create(chromosome, position)
                    targetRegionEnrichment.mTargetedRegions.add(genomePosition)
                    if ((relativeEnrichment != null) && !relativeEnrichment.isNaN())
                    {
                        targetRegionEnrichment.mTargetRelativeEnrichment[genomePosition] = relativeEnrichment
                    }
                }
            }
            return targetRegionEnrichment
        }
    }
}