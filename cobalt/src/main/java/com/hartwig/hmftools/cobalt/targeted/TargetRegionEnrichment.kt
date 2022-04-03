package com.hartwig.hmftools.cobalt.targeted

import com.hartwig.hmftools.common.genome.position.GenomePosition
import com.hartwig.hmftools.common.genome.position.GenomePositions
import com.hartwig.hmftools.common.utils.FileWriterUtils
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
                // skip first line which is header
                reader.readLine()
                var line: String? = reader.readLine()
                while (line != null)
                {
                    val values = line.split(DELIMITER)

                    val chromosome: String = values[0].intern()
                    val position: Int = (values[1].toDouble()).toInt()
                    val relativeEnrichment: Double? = values[2].toDoubleOrNull()
                    val genomePosition = GenomePositions.create(chromosome, position)
                    targetRegionEnrichment.mTargetedRegions.add(genomePosition)
                    if (relativeEnrichment != null && !relativeEnrichment.isNaN())
                    {
                        targetRegionEnrichment.mTargetRelativeEnrichment[genomePosition] = relativeEnrichment
                    }

                    line = reader.readLine()
                }
            }
            return targetRegionEnrichment
        }
    }
}