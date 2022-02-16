package com.hartwig.hmftools.teal.breakend

import kotlin.math.min

object BreakEndEvidenceFilter
{
    const val MAX_GERMLINE_SUPPORT_COUNT = 5
    const val MAX_GERMLINE_SUPPORT_RATIO = 0.02
    const val MIN_TELOMERIC_LENGTH = 20
    const val MIN_ANCHOR_LENGTH = 50
    const val MIN_SPLIT_READ_COUNT = 3
    const val MIN_DISCORDANT_PAIR_COUNT = 1
    const val MAX_COHORT_FREQ = 1000000
    const val MIN_TUMOR_MAPQ = 300

    private val supportSplitReadTypes = arrayOf(
        Pair(Fragment.AlignedReadType.SPLIT_READ_TELOMERIC, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC),
        Pair(Fragment.AlignedReadType.SPLIT_READ_TELOMERIC, Fragment.PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC),
        Pair(Fragment.AlignedReadType.SPLIT_READ_TELOMERIC, Fragment.PairedReadType.NOT_DISCORDANT),
        Pair(Fragment.AlignedReadType.SPLIT_READ_NOT_TELOMERIC, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC))

    private val supportDiscordantPairTypes = arrayOf(
        Pair(Fragment.AlignedReadType.DISCORDANT_PAIR, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC),
        Pair(Fragment.AlignedReadType.SPLIT_READ_TELOMERIC, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC),
        Pair(Fragment.AlignedReadType.SPLIT_READ_NOT_TELOMERIC, Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC))

    fun getFilterReasons(evidence: TelomericBreakEndEvidence) : List<String>
    {
        val filterList = ArrayList<String>()

        if (evidence.tumorSupport != null)
        {
            if (evidence.tumorSupport.key.isDuplicate)
            {
                filterList.add("duplicate")
            }

            if (evidence.tumorSupport.longestTelomereSegment == null ||
                evidence.tumorSupport.longestTelomereSegment!!.length < MIN_TELOMERIC_LENGTH)
            {
                filterList.add("minTelomericLength")
            }

            if (evidence.tumorSupport.longestSplitReadAlignLength < MIN_ANCHOR_LENGTH)
            {
                filterList.add("minAnchorLength")
            }

            if (evidence.tumorSupport.fragmentCount({ f -> supportSplitReadTypes.contains(Pair(f.alignedReadType, f.pairedReadType)) }) < MIN_SPLIT_READ_COUNT)
            {
                filterList.add("minSplitReads")
            }

            if (evidence.tumorSupport.fragmentCount({ f -> supportDiscordantPairTypes.contains(Pair(f.alignedReadType, f.pairedReadType)) }) < MIN_DISCORDANT_PAIR_COUNT)
            {
                filterList.add("minDiscordantPairs")
            }

            if (evidence.tumorSupport.sumMapQ({ f -> Pair(f.alignedReadType, f.pairedReadType) in Fragment.supportFragTypes }) < MIN_TUMOR_MAPQ)
            {
                filterList.add("minTumorMAPQ")
            }
        }
        else
        {
            filterList.add("notInTumor")
        }

        if (evidence.germlineSupport != null)
        {
            val tumorSupportFragCount = evidence.tumorSupport?.supportFragmentCount() ?: 0
            if (evidence.germlineSupport.supportFragmentCount() >
                min(MAX_GERMLINE_SUPPORT_COUNT.toDouble(), tumorSupportFragCount * MAX_GERMLINE_SUPPORT_RATIO))
            {
                filterList.add("maxGermlineSupport")
            }
        }

        if (filterList.isEmpty())
            filterList.add("PASS")

        return filterList
    }
}
