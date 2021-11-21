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
}

data class Fragment(
    val type: Type,
    val readGroup: ReadGroup,
    val splitReads: Collection<SAMRecord>,
    val facingBreakReads: Collection<SAMRecord>)
{
    enum class Type(val code: String)
    {
        SPLIT_READ("SR"),
        DISCORDANT_PAIR("DP"),
        CONTRA_SPLIT_READ("CSR"),
        CONTRA_DISCORDANT_PAIR("CDP")
    }
}

// an added telomere usually means a section of the chromosome is
// broken off
// there might be a case for working out if the telomere is attached to the centremere
// or somewhere else??
// there are two types of such locations
// 1. On the sen
data class TelomericBreakEndKey(
    // +1 orientation right side been lost and got a C rich telomere
    // -1 orientation left side been lost and got a G rich telomere
    val type: TelomericBreakEndType,
    val chromosome: String,
    var position: Int
) : Comparable<TelomericBreakEndKey>
{
    // Overriding compareTo() method to allow us to sort
    override fun compareTo(other: TelomericBreakEndKey): Int
    {
        return Comparator.comparing { obj: TelomericBreakEndKey -> obj.type }
            .thenComparing { obj: TelomericBreakEndKey -> obj.chromosome }
            .thenComparingInt { obj: TelomericBreakEndKey -> obj.position }
            .compare(this, other)
    }

    override fun toString(): String = "type(${type}) ${chromosome}:${position}"
}

data class TelomericBreakEnd(
    val key: TelomericBreakEndKey,
    val telomericSplitReads: MutableList<SAMRecord> = ArrayList(),
    val fragments: MutableList<Fragment> = ArrayList())
{
    constructor(type: TelomericBreakEndType,
                chromosome: String,
                position: Int
    ) : this(TelomericBreakEndKey(type, chromosome, position))

    val type: TelomericBreakEndType get() = key.type
    val chromosome: String get() = key.chromosome
    val position: Int get() = key.position

    fun fragmentCount(fragType: Fragment.Type) : Int
    {
        return fragments.count({ fragment -> fragment.type == fragType })
    }

    fun sumMapQ() : Int
    {
        var mapqSum = 0
        for (fragment in fragments)
        {
            mapqSum += if (fragment.splitReads.isNotEmpty())
            {
                fragment.splitReads.sumOf({ r -> r.mappingQuality })
            }
            else
            {
                fragment.facingBreakReads.sumOf({ r -> r.mappingQuality })
            }
        }
        return mapqSum
    }

    fun isRightTelomeric() : Boolean = type.isRightTelomeric()

    override fun toString(): String = "type(${type}) ${chromosome}:${position}"
}
