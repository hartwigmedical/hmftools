package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import java.io.File

object Cdr3Utils
{
    const val DEFAULT_PARTITION_SIZE = 10000000

    fun openSamReader(config: Cdr3Params): SamReader
    {
        var factory = SamReaderFactory.makeDefault()
        if (config.RefGenomePath != null && !config.RefGenomePath!!.isEmpty())
        {
            factory = factory.referenceSequence(File(config.RefGenomePath!!))
        }
        return factory.open(File(config.BamPath))
    }

    @JvmStatic
    fun createPartitions(config: Cdr3Params): List<GenomeRegion>
    {
        val samReader = openSamReader(config)
        val samSequences = samReader.fileHeader.sequenceDictionary.sequences
        val partitions: MutableList<GenomeRegion> = ArrayList()
        val partitionSize = DEFAULT_PARTITION_SIZE
        for (seq in samSequences)
        {
            val chrStr = seq.sequenceName
            if (config.SpecificChr != null && config.SpecificChr != chrStr)
                continue
            val chromosomeLength = seq.sequenceLength
            var startPos = 0
            while (startPos < chromosomeLength)
            {
                var endPos = startPos + partitionSize - 1
                if (endPos + partitionSize * 0.2 > chromosomeLength) endPos = chromosomeLength
                partitions.add(GenomeRegions.create(chrStr, startPos, endPos))
                startPos = endPos + 1
            }
        }
        return partitions
    }

    @JvmStatic
    fun countsToString(counts: IntArray): String
    {
        val RADIX = 36
        // create a support string
        val stringBuilder = StringBuilder(counts.size)

        for (count in counts)
        {
            if (count >= RADIX)
                stringBuilder.append('#')
            else
                stringBuilder.append(count.toString(radix = RADIX))
        }
        return stringBuilder.toString()
    }

    fun insertDashes(str: String, vararg dashPos: Int): String
    {
        val stringBuilder = StringBuilder(str)

        // we must go in reverse order
        for (i in dashPos.size - 1 downTo 0)
        {
            val pos = dashPos[i]
            if (pos < str.length && pos != 0)
            {
                stringBuilder.insert(pos, '-')
            }
        }

        return stringBuilder.toString()
    }

}