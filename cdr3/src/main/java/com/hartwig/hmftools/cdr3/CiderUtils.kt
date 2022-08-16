package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import java.io.File

object CiderUtils
{
    const val DEFAULT_PARTITION_SIZE = 10000000

    fun openSamReader(config: CiderParams): SamReader
    {
        var factory = SamReaderFactory.makeDefault()
        if (config.RefGenomePath != null && !config.RefGenomePath!!.isEmpty())
        {
            factory = factory.referenceSequence(File(config.RefGenomePath!!))
        }
        return factory.open(File(config.BamPath))
    }

    @JvmStatic
    fun createPartitions(config: CiderParams): List<GenomeRegion>
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

    fun conservedAA(vjGeneType: VJGeneType): Char
    {
        if (vjGeneType.vj == VJ.V)
            return 'C'
        if (vjGeneType == VJGeneType.IGHJ)
            return 'W'
        return 'F'
    }

    fun getPairedVjGeneType(vjGene: VJGeneType) : VJGeneType
    {
        return when (vjGene)
        {
            VJGeneType.IGHV -> VJGeneType.IGHJ
            VJGeneType.IGHJ -> VJGeneType.IGHV
            VJGeneType.IGKV -> VJGeneType.IGKJ
            VJGeneType.IGKJ -> VJGeneType.IGKV
            VJGeneType.IGLV -> VJGeneType.IGLJ
            VJGeneType.IGLJ -> VJGeneType.IGLV
            VJGeneType.TRAV -> VJGeneType.TRAJ
            VJGeneType.TRAJ -> VJGeneType.TRAV
            VJGeneType.TRBV -> VJGeneType.TRBJ
            VJGeneType.TRBJ -> VJGeneType.TRBV
            VJGeneType.TRDV -> VJGeneType.TRDJ
            VJGeneType.TRDJ -> VJGeneType.TRDV
            VJGeneType.TRGV -> VJGeneType.TRGJ
            VJGeneType.TRGJ -> VJGeneType.TRGV
        }
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

    // insert dashes. Note that this works even if some dash positions
    // are invalid
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

    fun safeSubstring(str: String, start: Int, end: Int): String
    {
        if (start >= str.length)
        {
            return ""
        }
        if (end >= str.length)
        {
            return str.substring(start)
        }
        return str.substring(start, end)
    }

    fun safeSubstring(str: String, range: IntRange): String
    {
        if (str.length <= range.first)
        {
            return ""
        }
        if (str.length <= range.last)
        {
            return str.substring(range.first)
        }
        return str.substring(range)
    }
}