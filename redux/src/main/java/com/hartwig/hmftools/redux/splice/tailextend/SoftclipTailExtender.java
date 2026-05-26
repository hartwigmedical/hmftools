package com.hartwig.hmftools.redux.splice.tailextend;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MISMATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SOFTCLIP;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.redux.splice.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.redux.splice.rescue.CigarShape;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;

// Recovers ref-matching bases that bwa-mem2 softclipped at the read tails. Runs after junction
// rescue so rescue's boundary lookups see bwa's original cigar; this pass cleans up what rescue
// didn't merge. Walks each terminal softclip against the ref, converts S->M while running
// mismatches stay below max(1, length / 10), and caps the extension at MaxExtension. Skips when
// an annotated intron starts within MaxExtension of the S boundary — that range belongs to the
// junction-rescue path.
public class SoftclipTailExtender
{
    private final RefSequenceSource mRefSource;
    private final AnnotatedJunctionIndex mJunctionGuard;
    private final TailExtensionConfig mConfig;
    private final TailExtensionStatistics mStatistics;

    public SoftclipTailExtender(
            final RefSequenceSource refSource, final AnnotatedJunctionIndex junctionGuard,
            final TailExtensionConfig config)
    {
        mRefSource = refSource;
        mJunctionGuard = junctionGuard;
        mConfig = config;
        mStatistics = new TailExtensionStatistics();
    }

    public TailExtensionStatistics statistics()
    {
        return mStatistics;
    }

    public TailExtensionResult tryExtend(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(!mConfig.Enabled || mRefSource == null || chromosome == null
                || cigar == null || readBases == null || readBases.length == 0)
            return TailExtensionResult.unchanged();

        final List<CigarShape.Element> elements = CigarShape.parse(cigar);
        if(elements.isEmpty())
            return TailExtensionResult.unchanged();
        if(CigarShape.hasHardClip(elements))
        {
            mStatistics.countSkippedComplexShape();
            return TailExtensionResult.unchanged();
        }

        mStatistics.countEvaluated();

        // trailing extension runs first because the leading-side math reads alignmentStart, not
        // the cigar tail; this keeps the two extensions independent on the same working list.
        final List<CigarShape.Element> working = new ArrayList<>(elements);

        final int trailExtended = tryExtendTrailing(chromosome, alignmentStart, readBases, working);
        final int leadExtended = tryExtendLeading(chromosome, alignmentStart, readBases, working);

        if(trailExtended == 0 && leadExtended == 0)
            return TailExtensionResult.unchanged();

        mStatistics.countExtended(leadExtended, trailExtended);
        final int newStart = alignmentStart - leadExtended;
        return new TailExtensionResult(true, newStart, CigarShape.format(working),
                leadExtended, trailExtended);
    }

    private int tryExtendTrailing(
            final String chromosome, final int alignmentStart, final byte[] readBases,
            final List<CigarShape.Element> working)
    {
        final int last = working.size() - 1;
        if(last < 1)
            return 0;

        final CigarShape.Element softclip = working.get(last);
        if(softclip.Op != OP_SOFTCLIP || softclip.Length < mConfig.MinSoftclipLength)
            return 0;

        final CigarShape.Element adjacent = working.get(last - 1);
        if(!isMatchedOp(adjacent.Op))
        {
            mStatistics.countSkippedComplexShape();
            return 0;
        }

        final int alignmentEnd = alignmentStart + CigarShape.referenceSpan(working) - 1;

        if(crossesAnnotatedJunctionTrailing(chromosome, alignmentEnd))
            return 0;

        final int extendBudget = Math.min(softclip.Length, mConfig.MaxExtension);
        final byte[] refBases = mRefSource.getBases(chromosome, alignmentEnd + 1, alignmentEnd + extendBudget);
        if(refBases == null || refBases.length < mConfig.MinExtension)
        {
            mStatistics.countSkippedNoRef();
            return 0;
        }

        final int walkLength = Math.min(refBases.length, extendBudget);
        final int readStart = readBases.length - softclip.Length;
        final int bestLength = longestAcceptedPrefix(readBases, readStart, refBases, 0, walkLength, true);

        if(bestLength < mConfig.MinExtension)
        {
            if(bestLength == 0)
                mStatistics.countRejectedTooManyMismatches();
            return 0;
        }

        applyExtension(working, last, last - 1, bestLength);
        return bestLength;
    }

    private int tryExtendLeading(
            final String chromosome, final int alignmentStart, final byte[] readBases,
            final List<CigarShape.Element> working)
    {
        if(working.size() < 2)
            return 0;

        final CigarShape.Element softclip = working.get(0);
        if(softclip.Op != OP_SOFTCLIP || softclip.Length < mConfig.MinSoftclipLength)
            return 0;

        final CigarShape.Element adjacent = working.get(1);
        if(!isMatchedOp(adjacent.Op))
        {
            mStatistics.countSkippedComplexShape();
            return 0;
        }

        if(crossesAnnotatedJunctionLeading(chromosome, alignmentStart))
            return 0;

        int extendBudget = Math.min(softclip.Length, mConfig.MaxExtension);
        if(alignmentStart - extendBudget < 1)
            extendBudget = alignmentStart - 1;
        if(extendBudget < mConfig.MinExtension)
            return 0;

        final int refStart = alignmentStart - extendBudget;
        final int refEnd = alignmentStart - 1;
        final byte[] refBases = mRefSource.getBases(chromosome, refStart, refEnd);
        if(refBases == null || refBases.length < mConfig.MinExtension)
        {
            mStatistics.countSkippedNoRef();
            return 0;
        }

        final int walkLength = Math.min(refBases.length, extendBudget);
        // walk backward from the M boundary into the softclip: read[softclipLen-1-i] vs
        // ref[refLen-1-i].
        final int bestLength = longestAcceptedPrefix(
                readBases, softclip.Length - 1, refBases, refBases.length - 1, walkLength, false);

        if(bestLength < mConfig.MinExtension)
        {
            if(bestLength == 0)
                mStatistics.countRejectedTooManyMismatches();
            return 0;
        }

        applyExtension(working, 0, 1, bestLength);
        return bestLength;
    }

    // longest p in [0, walkLength] with cumulative mismatches <= max(1, p / 10). The min-of-1
    // tolerance lets short tails (3-9 bp) recover a single mismatch.
    private static int longestAcceptedPrefix(
            final byte[] readBases, final int readAnchor,
            final byte[] refBases, final int refAnchor,
            final int walkLength, final boolean forward)
    {
        int mismatches = 0;
        int bestLength = 0;
        for(int i = 0; i < walkLength; ++i)
        {
            final int step = forward ? i : -i;
            final byte readBase = readBases[readAnchor + step];
            final byte refBase = refBases[refAnchor + step];
            if(!basesEqualIgnoreCase(readBase, refBase))
                ++mismatches;
            final int length = i + 1;
            final int allowed = Math.max(1, length / 10);
            if(mismatches <= allowed)
                bestLength = length;
            else if(mismatches > allowed + 1)
                break;
        }
        return bestLength;
    }

    private boolean crossesAnnotatedJunctionTrailing(final String chromosome, final int alignmentEnd)
    {
        if(mJunctionGuard == null)
            return false;
        for(int offset = 1; offset <= mConfig.MaxExtension; ++offset)
        {
            if(!mJunctionGuard.introByStart(chromosome, alignmentEnd + offset).isEmpty())
            {
                mStatistics.countSkippedForJunctionGuard();
                return true;
            }
        }
        return false;
    }

    private boolean crossesAnnotatedJunctionLeading(final String chromosome, final int alignmentStart)
    {
        if(mJunctionGuard == null)
            return false;
        for(int offset = 1; offset <= mConfig.MaxExtension; ++offset)
        {
            if(!mJunctionGuard.introByEnd(chromosome, alignmentStart - offset).isEmpty())
            {
                mStatistics.countSkippedForJunctionGuard();
                return true;
            }
        }
        return false;
    }

    private static boolean isMatchedOp(final char op)
    {
        return op == OP_MATCH || op == OP_SEQ_MATCH || op == OP_SEQ_MISMATCH;
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
            return true;
        return (a & ~0x20) == (b & ~0x20);
    }

    private static void applyExtension(
            final List<CigarShape.Element> working, final int softclipIndex,
            final int matchedIndex, final int extend)
    {
        final CigarShape.Element softclip = working.get(softclipIndex);
        final CigarShape.Element matched = working.get(matchedIndex);

        final CigarShape.Element grownMatch = new CigarShape.Element(matched.Length + extend, matched.Op);
        working.set(matchedIndex, grownMatch);

        if(softclip.Length - extend == 0)
        {
            working.remove(softclipIndex);
        }
        else
        {
            working.set(softclipIndex, new CigarShape.Element(softclip.Length - extend, OP_SOFTCLIP));
        }
    }
}
