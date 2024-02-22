package com.hartwig.hmftools.sage.evidence;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import htsjdk.samtools.SAMRecord;

public class ArtefactContext
{
    // if a variant is adjacent (factoring in left-alignment for INDELs) to an immediately upstream homopolymer of length >= 8,
    // cap per-read base qual to the lowest base qual from the first N bases of the homopolymer
    // Allow for right-alignment of DELs in the search for the start of the homopolymer

    private final int[] mHomopolymerStartOffset;
    private final boolean mRequiresCheck;

    public static final byte NOT_APPLICABLE_BASE_QUAL = -1;
    public static final byte NO_BASE = 0;

    private static final int NO_INDEX = -1;
    private static final int HOMOPOLYMER_BASE_SEARCH = 3;
    private static final int HOMOPOLYMER_REPEAT_LENGTH = 8;

    public ArtefactContext(final int[] homopolymerStartOffset)
    {
        mHomopolymerStartOffset = homopolymerStartOffset;
        mRequiresCheck = mHomopolymerStartOffset[SE_START] != NO_BASE || mHomopolymerStartOffset[SE_END] != NO_BASE;
    }

    public static ArtefactContext buildContext(final SimpleVariant variant, final IndexedBases indexedBases)
    {
        // check for any single-base repeat of 8+ bases not covering the variant
        String flankAndBases = indexedBases.fullString();

        int indexInBases = indexedBases.Index - indexedBases.LeftFlankIndex;

        // for SNVs and MNVs, restrict the search to within the variant bases
        int hpStartOffset = NO_INDEX;
        int hpEndOffset = NO_INDEX;

        int altLength = variant.alt().length();

        for(int varIndex = indexInBases; varIndex < indexInBases + altLength; ++varIndex)
        {
            if(hasHomopolymerRepeat(flankAndBases, varIndex, true))
            {
                hpStartOffset = indexInBases - varIndex;
            }

            if(hpEndOffset != NO_INDEX)
                break;

            if(hasHomopolymerRepeat(flankAndBases, varIndex, false))
            {
                hpEndOffset = varIndex - indexInBases;
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
                    hpEndOffset = delStartIndex - indexInBases;
                }
            }
        }

        if(hpStartOffset == NO_INDEX && hpEndOffset == NO_INDEX)
            return null;

        return new ArtefactContext(new int[] { hpStartOffset, hpEndOffset} );
    }

    public boolean requiresCheck() { return mRequiresCheck; }
    public int homopolymerOffset(int seIndex) { return mHomopolymerStartOffset[seIndex]; }
    public boolean hasHomopolymerOffset(int seIndex) { return mHomopolymerStartOffset[seIndex] != NO_INDEX; }

    public byte findApplicableBaseQual(final SAMRecord record, int varReadIndex)
    {
        int homopolSide = record.getReadNegativeStrandFlag() ? SE_END : SE_START;

        if(mHomopolymerStartOffset[homopolSide] == NO_INDEX)
            return NOT_APPLICABLE_BASE_QUAL;

        return findHomopolymerBaseQual(
                record, varReadIndex, mHomopolymerStartOffset[homopolSide], homopolSide == SE_START);
    }

    private byte findHomopolymerBaseQual(final SAMRecord record, int varReadIndex, int hpOffset, boolean searchDown)
    {
        int hpStartIndex = searchDown ? varReadIndex - hpOffset : varReadIndex + hpOffset;

        int minBaseQual = 0;

        for(int i = 0; i < HOMOPOLYMER_BASE_SEARCH; ++i)
        {
            int hpIndex = hpStartIndex + (searchDown ? -i : i);

            if(i == 0)
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
        return format("hpOffsets start(%d) end(%d)", mHomopolymerStartOffset[0], mHomopolymerStartOffset[1]);
    }
}
