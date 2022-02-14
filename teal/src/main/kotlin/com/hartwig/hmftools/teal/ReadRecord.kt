package com.hartwig.hmftools.teal

import htsjdk.samtools.Cigar
import htsjdk.samtools.SAMFlag
import htsjdk.samtools.SAMRecord

class ReadRecord(
    val Id: String, val Chromosome: String, val PosStart: Int, val PosEnd: Int, val ReadBases: String, baseQualityString: String,
    cigar: Cigar?, insertSize: Int, flags: Int, mateChromosome: String, matePosStart: Int)
{
    val BaseQualityString: String
    val Length: Int // of bases
    val Cigar: Cigar?
    private var mFlags: Int
    private val mMateChromosome: String
    private val mMatePosStart: Int
    private val mFragmentInsertSize: Int
    var suppAlignment: String?
    private var mMapQuality: Short
    private var mHasTeloContent: Boolean
    private var mCompleteGroup: Boolean
    fun hasTeloContent(): Boolean
    {
        return mHasTeloContent
    }

    fun setTeloContent(toggle: Boolean)
    {
        mHasTeloContent = toggle
    }

    fun completeGroup(): Boolean
    {
        return mCompleteGroup
    }

    fun markCompleteGroup()
    {
        mCompleteGroup = true
    }

    fun orientation(): Byte
    {
        // first in pair has orientation of +1 if not reversed, and vice versa for the second in the pair
        return if (isFirstOfPair) if (!isReadReversed) 1 else (-1.toByte()).toByte() else if (isReadReversed) (-1.toByte()).toByte() else 1
    }

    fun flags(): Int
    {
        return mFlags
    }

    val isReadReversed: Boolean
        get() = mFlags and SAMFlag.READ_REVERSE_STRAND.intValue() != 0
    val isFirstOfPair: Boolean
        get() = mFlags and SAMFlag.FIRST_OF_PAIR.intValue() != 0
    val isDuplicate: Boolean
        get() = mFlags and SAMFlag.DUPLICATE_READ.intValue() != 0
    val isUnmapped: Boolean
        get() = mFlags and SAMFlag.READ_UNMAPPED.intValue() != 0
    val isMateUnmapped: Boolean
        get() = mFlags and SAMFlag.MATE_UNMAPPED.intValue() != 0

    fun hasSuppAlignment(): Boolean
    {
        return suppAlignment != null
    }

    fun setMapQuality(mapQuality: Short)
    {
        mMapQuality = mapQuality
    }

    fun mapQuality(): Short
    {
        return mMapQuality
    }

    fun fragmentInsertSize(): Int
    {
        return mFragmentInsertSize
    }

    fun mateChromosome(): String
    {
        return mMateChromosome
    }

    fun mateStartPosition(): Int
    {
        return mMatePosStart
    }

    fun matches(other: ReadRecord): Boolean
    {
        return Id == other.Id && Cigar.toString() == other.Cigar.toString() && PosStart == other.PosStart && PosEnd == other.PosEnd
    }

    fun setFlag(flag: SAMFlag, toggle: Boolean)
    {
        mFlags = if (toggle) mFlags or flag.intValue() else mFlags and flag.intValue().inv()
    }

    override fun toString(): String
    {
        return String.format(
            "range(%s: %d -> %d) length(%d) cigar(%s) id(%s)",
            Chromosome, PosStart, PosEnd, Length, Cigar?.toString() ?: "", Id)
    }

    companion object
    {
        private const val SUPPLEMENTARY_ATTRIBUTE = "SA"
        fun from(record: SAMRecord): ReadRecord
        {
            val readId = if (record.isSecondaryAlignment) String.format(
                "%s_%s",
                record.readName, record.getAttribute("HI"))
            else record.readName
            val read = ReadRecord(
                readId, record.referenceName, record.start, record.end,
                record.readString, record.baseQualityString, record.cigar, record.inferredInsertSize, record.flags,
                record.mateReferenceName, record.mateAlignmentStart)
            read.suppAlignment = record.getStringAttribute(SUPPLEMENTARY_ATTRIBUTE)
            read.setMapQuality(record.mappingQuality.toShort())
            return read
        }
    }

    init
    {
        Length = ReadBases.length
        Cigar = cigar
        BaseQualityString = baseQualityString
        mFlags = flags
        mMateChromosome = mateChromosome
        mMatePosStart = matePosStart
        mFragmentInsertSize = insertSize
        suppAlignment = null
        mMapQuality = 0
        mHasTeloContent = false
        mCompleteGroup = false
    }
}