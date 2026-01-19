package com.hartwig.hmftools.teal

import com.hartwig.hmftools.teal.telbam.TelbamParams
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import htsjdk.samtools.ValidationStringency
import org.apache.commons.lang3.StringUtils
import java.io.File

object TealUtils
{
    fun hasTelomericContent(readBases: String): Boolean
    {
        for (teloSeq in TealConstants.CANONICAL_TELOMERE_SEQUENCES)
        {
            val matchIndex = readBases.indexOf(teloSeq)
            if (matchIndex != -1)
            {
                return true
            }
        }
        return false
    }

    // todo: try SequenceUtil.reverseComplement
    fun reverseComplementSequence(seq: String): String
    {
        val builder = StringBuilder()
        for (i in seq.length - 1 downTo 0)
        {
            val b = seq[i]
            when (b)
            {
                'A' ->
                {
                    builder.append('T')
                }
                'T' ->
                {
                    builder.append('A')
                }
                'G' ->
                {
                    builder.append('C')
                }
                'C' ->
                {
                    builder.append('G')
                }
            }
        }
        return builder.toString()
    }

    // determine if the read is poly G
    // note: We must supply this function with the read as is from
    // the read
    fun isPolyGC(readSeq: String): Boolean
    {
        // we use threshold of 90%
        val numGCs = Math.max(StringUtils.countMatches(readSeq, 'G'), StringUtils.countMatches(readSeq, 'C')).toDouble()
        val gcFrac = numGCs / readSeq.length
        return gcFrac >= TealConstants.POLY_G_THRESHOLD
    }

    fun openSamReader(config: TelbamParams): SamReader
    {
        var factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT)
        if (config.refGenomeFile != null && !config.refGenomeFile!!.isEmpty())
        {
            factory = factory.referenceSequence(File(config.refGenomeFile!!))
        }
        return factory.open(File(config.bamFile!!))
    }
}