package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.sage.quality.QualityCalculator.INVALID_BASE_QUAL;

import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class ArtefactContext
{
    // if a variant is adjacent (factoring in left-alignment for INDELs) to an immediately upstream homopolymer of length >= 8,
    // cap per-read base qual to the lowest base qual from the first N bases of the homopolymer
    // Allow for right-alignment of DELs in the search for the start of the homopolymer

    private final int[] mHomopolymerStartOffset;
    private final boolean mRequiresCheck;
    private final byte[] mHomopolymerBase;

    public static final byte NO_BASE = 0;

    private static final int NO_INDEX = -10; // beyond the permitted range
    private static final int HOMOPOLYMER_BASE_SEARCH = 3;
    private static final int HOMOPOLYMER_REPEAT_LENGTH = 8;

    public ArtefactContext(final int[] homopolymerStartOffset, final byte[] homopolymerBase)
    {
        mHomopolymerStartOffset = homopolymerStartOffset;
        mHomopolymerBase = homopolymerBase;
        mRequiresCheck = mHomopolymerStartOffset[SE_START] != NO_BASE || mHomopolymerStartOffset[SE_END] != NO_BASE;
    }

    public static ArtefactContext buildContext(final VariantReadContext readContext)
    {
        // check for any single-base repeat of 8+ bases not covering the variant
        String flankAndBases = readContext.readBases();

        int indexInBases = readContext.VarIndex;

        SimpleVariant variant = readContext.variant();

        // for SNVs and MNVs, restrict the search to within the variant bases
        int hpStartOffset = NO_INDEX;
        int hpEndOffset = NO_INDEX;
        byte hpStartBase = NO_INDEX;
        byte hpEndBase = NO_INDEX;

        int altLength = variant.alt().length();

        for(int varIndex = indexInBases; varIndex < indexInBases + altLength; ++varIndex)
        {
            if(hasHomopolymerRepeat(flankAndBases, varIndex, true))
            {
                hpStartBase = (byte)flankAndBases.charAt(varIndex);
                hpStartOffset = indexInBases - varIndex;

                // check for a single homopolymer base downstream of the variant eg from an insert
                if(flankAndBases.charAt(varIndex + 1) == flankAndBases.charAt(varIndex))
                    --hpStartOffset;
            }

            if(hpEndOffset != NO_INDEX)
                break;

            if(hasHomopolymerRepeat(flankAndBases, varIndex, false))
            {
                hpEndBase = (byte)flankAndBases.charAt(varIndex);
                hpEndOffset = varIndex - indexInBases;

                if(flankAndBases.charAt(varIndex - 1) == flankAndBases.charAt(varIndex))
                    --hpEndOffset;
            }
            else if(variant.isDelete())
            {
                String deletedBases = variant.ref().substring(1);
                int delBaseLength = deletedBases.length();
                int delStartIndex = indexInBases + 1;

                while(delStartIndex + delBaseLength < flankAndBases.length()
                && deletedBases.equals(flankAndBases.substring(delStartIndex, delStartIndex + delBaseLength)))
                {
                    delStartIndex += deletedBases.length();
                }

                if(hasHomopolymerRepeat(flankAndBases, delStartIndex, false))
                {
                    hpEndBase = (byte)flankAndBases.charAt(delStartIndex);
                    hpEndOffset = delStartIndex - indexInBases;
                }
            }
        }

        if(hpStartOffset == NO_INDEX && hpEndOffset == NO_INDEX)
            return null;

        return new ArtefactContext(new int[] { hpStartOffset, hpEndOffset}, new byte[] { hpStartBase, hpEndBase });
    }

    public boolean requiresCheck() { return mRequiresCheck; }
    public int homopolymerOffset(int seIndex) { return mHomopolymerStartOffset[seIndex]; }
    public int homopolymerBase(int seIndex) { return mHomopolymerBase[seIndex]; }
    public boolean hasHomopolymerOffset(int seIndex) { return mHomopolymerStartOffset[seIndex] != NO_INDEX; }

    public byte findApplicableBaseQual(final SAMRecord record, int varReadIndex)
    {
        int homopolymerSide = record.getReadNegativeStrandFlag() ? SE_END : SE_START;

        if(mHomopolymerStartOffset[homopolymerSide] == NO_INDEX)
            return INVALID_BASE_QUAL;

        return findHomopolymerBaseQual(
                record, varReadIndex, mHomopolymerStartOffset[homopolymerSide], mHomopolymerBase[homopolymerSide],
                homopolymerSide == SE_START);
    }

    private byte findHomopolymerBaseQual(final SAMRecord record, int varReadIndex, int hpOffset, byte hpBase, boolean searchDown)
    {
        int hpStartIndex = searchDown ? varReadIndex - hpOffset : varReadIndex + hpOffset;

        int minBaseQual = -1;

        for(int i = 0; i < HOMOPOLYMER_BASE_SEARCH; ++i)
        {
            int hpIndex = hpStartIndex + (searchDown ? -i : i);

            if(hpIndex < 0 || hpIndex >= record.getBaseQualities().length)
                break;

            if(record.getReadBases()[hpIndex] != hpBase)
                continue;

            if(minBaseQual < 0)
                minBaseQual = record.getBaseQualities()[hpIndex];
            else
                minBaseQual = min(minBaseQual, record.getBaseQualities()[hpIndex]);
        }

        return (byte)minBaseQual;
    }

    private static boolean hasHomopolymerRepeat(final String bases, int startIndex, boolean searchDown)
    {
        int repeatCount = 1;
        char lastChar = bases.charAt(startIndex);

        int i = startIndex + (searchDown ? -1 : 1);

        while(i >= 0 && i < bases.length())
        {
            char nextChar = bases.charAt(i);
            if(lastChar == nextChar)
            {
                ++repeatCount;

                if(repeatCount >= HOMOPOLYMER_REPEAT_LENGTH)
                    return true;
            }
            else
            {
                return false;
            }

            if(searchDown)
                --i;
            else
                ++i;
        }

        return false;
    }

    public String toString()
    {
        return format("bases(%c / %c) offsets(%d / %d)",
                (char)mHomopolymerBase[0], (char)mHomopolymerBase[1], mHomopolymerStartOffset[0], mHomopolymerStartOffset[1]);
    }
}
