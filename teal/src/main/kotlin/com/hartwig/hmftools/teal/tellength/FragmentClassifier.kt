package com.hartwig.hmftools.teal.tellength

import com.hartwig.hmftools.common.sequencing.SequencingType
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils
import org.apache.commons.lang3.StringUtils
import org.apache.logging.log4j.LogManager
import java.util.regex.Matcher
import java.util.regex.Pattern

class FragmentClassifier(val sequencingType: SequencingType)
{
    private val minCanonicalCount: Int
    private val minConsecutiveHexamers: Int
    private val singleReadTelomereCheckBases: Int

    // make matcher thread local since they are not thread safe
    private val gTeloPatternMatcher: ThreadLocal<Matcher>
    private val cTeloPatternMatcher: ThreadLocal<Matcher>

    init
    {
        when(sequencingType) {
            SequencingType.ILLUMINA -> {
                minCanonicalCount = 4
                minConsecutiveHexamers = 6
                singleReadTelomereCheckBases = 151
            }
            SequencingType.ULTIMA -> {
                minCanonicalCount = 3
                minConsecutiveHexamers = 0
                singleReadTelomereCheckBases = 75
            }
            SequencingType.SBX -> {
                minCanonicalCount = 3
                minConsecutiveHexamers = 0
                singleReadTelomereCheckBases = 75
            }
        }

        // we match the start of the sequence, must be at least 6 repeats of telomere. The reason we only match the start of the sequence
        // is to account for poly G tail that many reads have
        val gTeloPattern =
            Pattern.compile("^[ACGT]{0,2}G{0,3}" + StringUtils.repeat("(?!GGGGGG)([ACGT][ACGT][ACGT]GGG)", minConsecutiveHexamers))
        val cTeloPattern =
            Pattern.compile("^[ACGT]{0,2}C{0,3}" + StringUtils.repeat("(?!CCCCCC)([ACGT][ACGT][ACGT]CCC)", minConsecutiveHexamers))

        gTeloPatternMatcher = ThreadLocal.withInitial { gTeloPattern.matcher("") }
        cTeloPatternMatcher = ThreadLocal.withInitial { cTeloPattern.matcher("") }
    }

    fun classifyFragment(readGroup: ReadGroup): FragmentType
    {
        var seq1: String
        var seq2: String

        if (readGroup.isPairedRead)
        {
            if (readGroup.reads.size < 2)
            {
                return FragmentType.UNKNOWN
            }

            // we want to check
            val read1 = readGroup.firstOfPair ?: return FragmentType.UNKNOWN
            val read2 = readGroup.secondOfPair ?: return FragmentType.UNKNOWN

            seq1 = read1.readString!!
            seq2 = read2.readString!!

            // qual_str1 = read1['BaseQualities']
            // qual_str2 = read2['BaseQualities']
            if (read1.readNegativeStrandFlag)
            {
                seq1 = TealUtils.reverseComplementSequence(seq1)
            }
            if (read2.readNegativeStrandFlag) {
                seq2 = TealUtils.reverseComplementSequence(seq2)
            }
        }
        else
        {
            // for single-end reads, we look at 75 bases from both ends. Take the reverse complement from the end to make it synonymous
            // with the paired-end case.

            val read = readGroup.allReads.first()
            val readBases = if (read.readNegativeStrandFlag) TealUtils.reverseComplementSequence(read.readString) else read.readString

            seq1 = readBases.take(singleReadTelomereCheckBases)
            seq2 = TealUtils.reverseComplementSequence(readBases).take(singleReadTelomereCheckBases)
        }

        // try each one in tern
        var fragType = classifyFragment(seq1, seq2)
        if (fragType == FragmentType.NOT_TELOMERE)
        {
            fragType = classifyFragment(seq2, seq1)
        }

        return fragType
    }

    private fun classifyFragment(readPairG: String, readPairC: String): FragmentType
    {
        val isGTelomeric = isLikelyGTelomeric(readPairG)
        val isCTelomeric = isLikelyCTelomeric(readPairC)

        // first we want to find the TTAGGG motif, then fill it up backwards
        if (isGTelomeric && isCTelomeric)
        {
            return FragmentType.F1
        }
        if (isGTelomeric)
        {
            // only g is telomeric
            return FragmentType.F4
        }
        return if (isCTelomeric)
        {
            // only C term is telomeric, could be F2a or F2b
            FragmentType.F2
        }
        else
        {
            // not telomeric
            FragmentType.NOT_TELOMERE
        }
    }

    fun isLikelyGTelomeric(readBases: String): Boolean
    {
        if (StringUtils.countMatches(readBases, "TTAGGG") >= minCanonicalCount)
        {
            val m = gTeloPatternMatcher.get()
            m.reset(readBases)
            return m.find()
        }
        return false
    }

    fun isLikelyCTelomeric(readBases: String): Boolean
    {
        if (StringUtils.countMatches(readBases, "CCCTAA") >= minCanonicalCount)
        {
            val m = cTeloPatternMatcher.get()
            m.reset(readBases)
            return m.find()
        }
        return false
    }

    companion object
    {
        private val logger = LogManager.getLogger(javaClass)
    }
}