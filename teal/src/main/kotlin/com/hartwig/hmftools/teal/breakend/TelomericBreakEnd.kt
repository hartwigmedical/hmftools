package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.teal.ReadGroup
import htsjdk.samtools.SAMRecord
import java.util.Comparator

enum class TelomericBreakEndType
{
    // potential new telomere
    RIGHT_G_TELOMERIC,
    LEFT_C_TELOMERIC,
    // potential join with telomere
    RIGHT_C_TELOMERIC,
    LEFT_G_TELOMERIC;

    fun isRightTelomeric() : Boolean
    {
        return this == RIGHT_G_TELOMERIC || this == RIGHT_C_TELOMERIC
    }

    fun isGTelomeric() : Boolean
    {
        return this == LEFT_G_TELOMERIC || this == RIGHT_G_TELOMERIC
    }
}

data class Fragment(
    val readGroup: ReadGroup,
    val alignedReadType: AlignedReadType,
    val pairedReadType: PairedReadType,
    val alignedRead: SAMRecord,
    val pairedRead: SAMRecord)
{
    enum class AlignedReadType(val code: String)
    {
        // split reads with split section telomeric
        SPLIT_READ_TELOMERIC("ALIGNED_SR"),
        // reads that span the breakend with no split
        NOT_SPLIT_READ("ALIGNED_NSR"),
        // split read where the split section is non telomeric
        SPLIT_READ_NOT_TELOMERIC("ALIGNED_SR_NTEL"),
        DISCORDANT_PAIR("ALIGNED_DP");

        fun isSplitRead() : Boolean
        {
            return this == SPLIT_READ_TELOMERIC || this == SPLIT_READ_NOT_TELOMERIC
        }
    }

    enum class PairedReadType(val code: String)
    {
        // split reads with split section telomeric
        DISCORDANT_PAIR_TELOMERIC("PAIRED_DP"),
        // discordant pair where the read is not telomeric
        DISCORDANT_PAIR_NOT_TELOMERIC("PAIRED_DP_NTEL"),
        // the pair is on the other side
        NOT_DISCORDANT("PAIRED_ND")
    }
    
    companion object
    {
        val supportFragTypes = arrayOf(
            Pair(AlignedReadType.SPLIT_READ_TELOMERIC, PairedReadType.DISCORDANT_PAIR_TELOMERIC),
            Pair(AlignedReadType.SPLIT_READ_TELOMERIC, PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC),
            Pair(AlignedReadType.SPLIT_READ_TELOMERIC, PairedReadType.NOT_DISCORDANT),
            Pair(AlignedReadType.SPLIT_READ_NOT_TELOMERIC, PairedReadType.DISCORDANT_PAIR_TELOMERIC),
            Pair(AlignedReadType.DISCORDANT_PAIR, PairedReadType.DISCORDANT_PAIR_TELOMERIC))

        val contraNotTelomericFragTypes = arrayOf(
            Pair(AlignedReadType.SPLIT_READ_NOT_TELOMERIC, PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC),
            Pair(AlignedReadType.SPLIT_READ_NOT_TELOMERIC, PairedReadType.NOT_DISCORDANT),
            Pair(AlignedReadType.DISCORDANT_PAIR, PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC),
            Pair(AlignedReadType.DISCORDANT_PAIR, PairedReadType.NOT_DISCORDANT))
    }
}

// an added telomere usually means a section of the chromosome is
// broken off
// there might be a case for working out if the telomere is attached to the centremere
// or somewhere else??
// there are two types of such locations
// 1. On the sen
data class TelomericBreakEnd(
    // +1 orientation right side been lost and got a C rich telomere
    // -1 orientation left side been lost and got a G rich telomere
    val type: TelomericBreakEndType,
    val chromosome: String,
    var position: Int,
    var isDuplicate: Boolean = false
) : Comparable<TelomericBreakEnd>
{
    // Overriding compareTo() method to allow us to sort
    override fun compareTo(other: TelomericBreakEnd): Int
    {
        return Comparator.comparing { obj: TelomericBreakEnd -> obj.type }
            .thenComparing { obj: TelomericBreakEnd -> obj.chromosome }
            .thenComparingInt { obj: TelomericBreakEnd -> obj.position }
            .compare(this, other)
    }

    override fun toString(): String = "type(${type}) ${chromosome}:${position} ${if (isDuplicate) "dup" else ""})"

    fun isRightTelomeric() : Boolean = type.isRightTelomeric()
    fun isGTelomeric() : Boolean = type.isGTelomeric()
}

data class TelomericBreakEndSupport(
    val key: TelomericBreakEnd,
    val fragments: MutableList<Fragment> = ArrayList(),
    var longestTelomereSegment: String? = null,
    var distanceToTelomere: Int? = null,
    var longestSplitReadAlignLength: Int = 0)
{
    constructor(type: TelomericBreakEndType,
                chromosome: String,
                position: Int
    ) : this(TelomericBreakEnd(type, chromosome, position))

    val type: TelomericBreakEndType get() = key.type
    val chromosome: String get() = key.chromosome
    val position: Int get() = key.position

    override fun toString(): String = "type(${type}) ${chromosome}:${position}"

    fun fragmentCount(fragmentFilter: (Fragment) -> Boolean) : Int
    {
        return fragments.count(fragmentFilter)
    }

    fun fragmentCount(alignedReadType: Fragment.AlignedReadType, pairedReadType: Fragment.PairedReadType) : Int
    {
        return fragments.count({ f -> f.alignedReadType == alignedReadType && f.pairedReadType == pairedReadType })
    }

    fun supportFragmentCount() : Int
    {
        return fragmentCount({ f -> Pair(f.alignedReadType, f.pairedReadType) in Fragment.supportFragTypes })
    }

    fun numSplitReads() : Int
    {
        return fragmentCount({ f -> f.alignedReadType.isSplitRead() })
    }

    fun sumMapQ(fragmentFilter: (Fragment) -> Boolean) : Int
    {
        var mapqSum = 0
        for (fragment in fragments)
        {
            if (fragmentFilter(fragment))
                mapqSum += fragment.alignedRead.mappingQuality
        }
        return mapqSum
    }
}

data class TelomericBreakEndEvidence(
    val tumorSupport: TelomericBreakEndSupport? = null,
    val germlineSupport: TelomericBreakEndSupport? = null)
{
    fun mainBreakEnd() : TelomericBreakEndSupport?
    {
        return tumorSupport ?: germlineSupport
    }

    fun chromosome() : String?
    {
        return mainBreakEnd()?.chromosome
    }

    fun position() : Int
    {
        return mainBreakEnd()?.position ?: 0
    }
}
