package com.hartwig.hmftools.teal

import java.lang.StringBuilder
import java.lang.ThreadLocal
import com.hartwig.hmftools.teal.telbam.TelbamParams
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import org.apache.commons.lang3.StringUtils
import java.io.File
import java.util.regex.Pattern

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

    private const val MIN_CANONICAL_COUNT = 4
    private const val MIN_CONSECUTIVE_HEXAMERS = 6

    // we match the start of the sequence, must be at least 6 repeats of telomere. The reason we only match the start of the sequence
    // is to account for poly G tail that many reads have
    private val sGTeloPattern =
        Pattern.compile("^[ACGT]{0,2}G{0,3}" + StringUtils.repeat("(?!GGGGGG)([ACGT][ACGT][ACGT]GGG)", MIN_CONSECUTIVE_HEXAMERS))
    private val sCTeloPattern =
        Pattern.compile("^[ACGT]{0,2}C{0,3}" + StringUtils.repeat("(?!CCCCCC)([ACGT][ACGT][ACGT]CCC)", MIN_CONSECUTIVE_HEXAMERS))

    // make matcher thread local since they are not thread safe
    private val sGTeloPatternMatcher = ThreadLocal.withInitial { sGTeloPattern.matcher("") }
    private val sCTeloPatternMatcher = ThreadLocal.withInitial { sCTeloPattern.matcher("") }
    fun isLikelyGTelomeric(readBases: String): Boolean
    {
        if (StringUtils.countMatches(readBases, "TTAGGG") >= MIN_CANONICAL_COUNT)
        {
            val m = sGTeloPatternMatcher.get()
            m.reset(readBases)
            return m.find()
        }
        return false
    }

    fun isLikelyCTelomeric(readBases: String): Boolean
    {
        if (StringUtils.countMatches(readBases, "CCCTAA") >= MIN_CANONICAL_COUNT)
        {
            val m = sCTeloPatternMatcher.get()
            m.reset(readBases)
            return m.find()
        }
        return false
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
        var factory = SamReaderFactory.makeDefault()
        if (config.refGenomeFile != null && !config.refGenomeFile!!.isEmpty())
        {
            factory = factory.referenceSequence(File(config.refGenomeFile!!))
        }
        return factory.open(File(config.bamFile))
    }
}