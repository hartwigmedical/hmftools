package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import com.hartwig.hmftools.sage.common.IndexedBases;

import htsjdk.samtools.SAMRecord;

public class ArtefactContext
{
    /*
    if a variant is adjacent (factoring in left-alignment for INDELs) to an immediately upstream homopolymer of length >= 8,
    cap per-read base qual to homopolymerRefBaseQual

    if the variant extends the homopolymer (specifically, deletes the first base downstream of the homopolymer, or an SNV with the homopolymer base as ALT),
    let N be the number of ref bases downstream of the deleted/SNV'd base(s) before a non-homopolymer base is encountered. then homopolymerRefBaseQual is
    the BQ of the (N+1)th last base of the homopolymer

    if the above applies, and there is a variant affecting the first homopolymer immediately upstream of the 8+ len homopolymer, or a variant extending
    or contracting the long homopolymer itself, cap its qual at the same homopolymerRefBaseQual
    */

    private final int mVariantSpan;
    private final byte[] mHomopolymerBases;
    private final int[] mSkipCounts; // 'N' in the notes above
    private final boolean mRequiresCheck;

    public static final byte NOT_APPLICABLE_BASE_QUAL = -1;
    public static final byte NO_BASE = 0;

    private static final int HOMOPOLYMER_REPEAT_LENGTH = 8;

    public ArtefactContext(final int variantSpan, final byte[] homopolymerBases, final int[] skipCounts)
    {
        mVariantSpan = variantSpan;
        mHomopolymerBases = homopolymerBases;
        mSkipCounts = skipCounts;
        mRequiresCheck = mHomopolymerBases[SE_START] != NO_BASE || mHomopolymerBases[SE_END] != NO_BASE;
    }

    public static ArtefactContext buildContext(final ReadContextCounter variant)
    {
        // check for any single-base repeat of 8+ bases not covering the variant
        final IndexedBases indexedBases = variant.readContext().indexedBases();
        String flankAndBases = indexedBases.fullString();

        int indexInBases = indexedBases.Index - indexedBases.LeftFlankIndex;
        int varIndexStart = variant.isIndel() ? indexInBases : indexInBases - 1; // move 1 down for an SNV/MNV
        int altLength = variant.variant().alt().length();
        int varIndexEnd = indexInBases + altLength;

        byte homopolymerBaseOnLeft = findHomopolymerBase(flankAndBases, varIndexStart, true);
        byte homopolymerBaseOnRight = findHomopolymerBase(flankAndBases, varIndexEnd, false);

        int rightAlignmentCount = 0;
        if(homopolymerBaseOnRight == NO_BASE && variant.variant().isDelete())
        {
            String deletedBases = variant.ref().substring(1);
            int delBaseLength = deletedBases.length();
            int delStartIndex = indexInBases + 1;
            while(deletedBases.equals(flankAndBases.substring(delStartIndex, delStartIndex + delBaseLength)))
            {
                delStartIndex += deletedBases.length();
                rightAlignmentCount += deletedBases.length();
            }

            homopolymerBaseOnRight = findHomopolymerBase(flankAndBases, delStartIndex, false);
        }

        if(homopolymerBaseOnLeft == NO_BASE && homopolymerBaseOnRight == NO_BASE)
            return null;

        int skipCountLeft = 0;
        int skipCountRight = 0;
        if(homopolymerBaseOnLeft != NO_BASE)
        {
            skipCountLeft = findSkipCount(flankAndBases, varIndexEnd, homopolymerBaseOnLeft, false);
        }

        if(homopolymerBaseOnRight != NO_BASE)
        {
            /*
            if(useRightAlignment)
                skipCountRight = 0;
            else
            */
            skipCountRight = findSkipCount(flankAndBases, varIndexStart + rightAlignmentCount, homopolymerBaseOnRight, true);
        }

        // calculate bases to get from the start of the variant to the first upstream ref base
        int variantSpan = variant.variant().isDelete() ? 2 : altLength;

        return new ArtefactContext(
                variantSpan, new byte[] { homopolymerBaseOnLeft, homopolymerBaseOnRight}, new int[] { skipCountLeft, skipCountRight} );
    }

    public boolean requiresCheck() { return mRequiresCheck; }
    public byte[] homopolymerBases() { return mHomopolymerBases; }
    public int[] skipCounts() { return mSkipCounts; }

    private static int findSkipCount(final String bases, final int startIndex, final byte homopolymerBase, boolean searchDown)
    {
        int skipCount = 0;

        int i = startIndex;

        while(i >= 0 && i < bases.length())
        {
            if((byte)bases.charAt(i) == homopolymerBase)
                ++skipCount;
            else
                break;

            if(searchDown)
                --i;
            else
                ++i;
        }

        return skipCount;
    }

    public byte findApplicableBaseQual(final ReadContextCounter variant, final SAMRecord record, int varReadIndex)
    {
        int homopolSide = record.getReadNegativeStrandFlag() ? SE_END : SE_START;

        if(mHomopolymerBases[homopolSide] == NO_BASE)
            return NOT_APPLICABLE_BASE_QUAL;

        return findHomopolymerBaseQual(
                record, varReadIndex, mHomopolymerBases[homopolSide], mSkipCounts[homopolSide], homopolSide == SE_END);
    }

    private byte findHomopolymerBaseQual(final SAMRecord record, int varReadIndex, byte homopolymerBase, int skipCount, boolean searchDown)
    {
        int baseIndex = searchDown ? varReadIndex + mVariantSpan : varReadIndex - 1;

        while(baseIndex >= 0 && baseIndex < record.getReadBases().length)
        {
            byte base = record.getReadBases()[baseIndex];

            if(base != homopolymerBase)
                break;

            if(searchDown)
                --baseIndex;
            else
                ++baseIndex;
        }

        int nBase = searchDown ? skipCount + 1 : -(skipCount + 1);

        int homopolymerBaseIndex = baseIndex + nBase;
        return record.getBaseQualities()[homopolymerBaseIndex];
    }

    private static byte findHomopolymerBase(final String bases, int startIndex, boolean searchDown)
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
                    return (byte)lastChar;
            }
            else
            {
                return NO_BASE;
            }

            if(searchDown)
                --i;
            else
                ++i;
        }

        return NO_BASE;
    }

    public String toString()
    {
        return format("bases(%d - %d) skip(%d - %d)", mHomopolymerBases[0], mHomopolymerBases[1], mSkipCounts[0], mSkipCounts[1]);
    }

    /*
    private static boolean withinLongHomoploymerRange(final ReadContextCounter readContextCounter)
    {
        // check for any single-base repeat of 8+ bases not covering the variant
        final IndexedBases indexedBases = readContextCounter.readContext().indexedBases();
        String flankAndBases = indexedBases.fullString();

        int varIndexStart = indexedBases.Index - indexedBases.LeftFlankIndex;

        int varIndexEnd = readContextCounter.variant().isDelete() ?
                varIndexStart + 1 : varIndexStart + readContextCounter.variant().alt().length() - 1;

        if(hasLongRepeat(flankAndBases, 0, varIndexStart - 1))
            return true;

        if(hasLongRepeat(flankAndBases, varIndexEnd + 1, flankAndBases.length() - 1))
            return true;

        return false;
    }
    */
}
