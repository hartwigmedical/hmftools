package com.hartwig.hmftools.tars.liftback.tailextend;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.tars.liftback.rescue.AnnotatedJunctionIndex;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Recovers ref-matching bases bwa-mem2 softclipped at read tails. Runs after junction rescue so
// rescue's boundary lookups see the original cigar. Skips when an annotated intron starts within
// MaxExtension of the softclip boundary - that range belongs to junction rescue.
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

        final List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.isEmpty())
            return TailExtensionResult.unchanged();
        if(CigarUtils.hasHardClip(elements))
        {
            mStatistics.countSkippedComplexShape();
            return TailExtensionResult.unchanged();
        }

        mStatistics.countEvaluated();

        final List<CigarElement> working = new ArrayList<>(elements);
        final int alignmentEnd = alignmentStart + CigarUtils.cigarAlignedLength(working) - 1;

        final int trailExtended = tryExtendSide(Side.TRAILING, chromosome, alignmentStart, alignmentEnd, readBases, working);
        final int leadExtended = tryExtendSide(Side.LEADING, chromosome, alignmentStart, alignmentEnd, readBases, working);

        if(trailExtended == 0 && leadExtended == 0)
            return TailExtensionResult.unchanged();

        mStatistics.countExtended(leadExtended, trailExtended);
        final int newStart = alignmentStart - leadExtended;
        return new TailExtensionResult(true, newStart, CigarUtils.cigarElementsToStr(working),
                leadExtended, trailExtended);
    }

    private int tryExtendSide(
            final Side side, final String chromosome, final int alignmentStart, final int alignmentEnd,
            final byte[] readBases, final List<CigarElement> working)
    {
        if(working.size() < 2)
            return 0;

        final int softclipIndex = side.softclipIndex(working);
        final int matchedIndex = side.matchedIndex(working);

        final CigarElement softclip = working.get(softclipIndex);
        if(softclip.getOperator() != CigarOperator.S || softclip.getLength() < mConfig.MinSoftclipLength)
            return 0;

        if(!isMatchedOp(working.get(matchedIndex).getOperator()))
        {
            mStatistics.countSkippedComplexShape();
            return 0;
        }

        final int boundary = side.refBoundary(alignmentStart, alignmentEnd);

        int extendBudget = Math.min(softclip.getLength(), mConfig.MaxExtension);
        if(side == Side.LEADING && alignmentStart - extendBudget < 1)
            extendBudget = alignmentStart - 1;
        if(extendBudget < mConfig.MinExtension)
            return 0;

        final byte[] refBases = side.fetchRef(mRefSource, chromosome, boundary, extendBudget);
        if(refBases == null || refBases.length < mConfig.MinExtension)
        {
            mStatistics.countSkippedNoRef();
            return 0;
        }

        final int walkLength = Math.min(refBases.length, extendBudget);
        final int bestLength = side.walk(readBases, refBases, softclip.getLength(), walkLength);

        // Junction guard applied after the walk: if the read matches the full window (intron retention),
        // extension is correct. Only defer when the walk stalls short, meaning the read diverges at the
        // donor and genuinely needs splicing.
        if(nearAnnotatedJunction(side, chromosome, boundary) && bestLength < walkLength)
        {
            mStatistics.countSkippedForJunctionGuard();
            return 0;
        }

        if(bestLength < mConfig.MinExtension)
        {
            if(bestLength == 0)
                mStatistics.countRejectedTooManyMismatches();
            return 0;
        }

        applyExtension(working, softclipIndex, matchedIndex, bestLength);
        return bestLength;
    }

    private boolean nearAnnotatedJunction(final Side side, final String chromosome, final int boundary)
    {
        if(mJunctionGuard == null)
            return false;
        for(int offset = 1; offset <= mConfig.MaxExtension; ++offset)
        {
            if(!side.junctionLookup(mJunctionGuard, chromosome, boundary, offset).isEmpty())
                return true;
        }
        return false;
    }

    private static boolean isMatchedOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    private static void applyExtension(
            final List<CigarElement> working, final int softclipIndex,
            final int matchedIndex, final int extend)
    {
        final CigarElement softclip = working.get(softclipIndex);
        final CigarElement matched = working.get(matchedIndex);

        final CigarElement grownMatch = new CigarElement(matched.getLength() + extend, matched.getOperator());
        working.set(matchedIndex, grownMatch);

        if(softclip.getLength() - extend == 0)
        {
            working.remove(softclipIndex);
        }
        else
        {
            working.set(softclipIndex, new CigarElement(softclip.getLength() - extend, CigarOperator.S));
        }
    }

    // Abstracts leading/trailing geometry so tryExtendSide is direction-free.
    private enum Side
    {
        LEADING
                {
                    @Override
                    int softclipIndex(final List<CigarElement> cigar) { return 0; }

                    @Override
                    int matchedIndex(final List<CigarElement> cigar) { return 1; }

                    @Override
                    int refBoundary(final int alignmentStart, final int alignmentEnd) { return alignmentStart; }

                    @Override
                    byte[] fetchRef(final RefSequenceSource src, final String chr, final int boundary, final int budget)
                    {
                        return src.getBases(chr, boundary - budget, boundary - 1);
                    }

                    @Override
                    int walk(final byte[] readBases, final byte[] refBases, final int softclipLength, final int walkLength)
                    {
                        // leading softclip's M boundary is at its inner end, so reverse to walk boundary-outward
                        final byte[] readWin = BoundaryReclaim.reversed(
                                Arrays.copyOfRange(readBases, softclipLength - walkLength, softclipLength));
                        final byte[] refWin = BoundaryReclaim.reversed(
                                Arrays.copyOfRange(refBases, refBases.length - walkLength, refBases.length));
                        return BoundaryReclaim.maxScoringPrefix(readWin, refWin);
                    }

                    @Override
                    List<ChrBaseRegion> junctionLookup(
                            final AnnotatedJunctionIndex idx, final String chr, final int boundary, final int offset)
                    {
                        return idx.introByEnd(chr, boundary - offset);
                    }
                },
        TRAILING
                {
                    @Override
                    int softclipIndex(final List<CigarElement> cigar) { return cigar.size() - 1; }

                    @Override
                    int matchedIndex(final List<CigarElement> cigar) { return cigar.size() - 2; }

                    @Override
                    int refBoundary(final int alignmentStart, final int alignmentEnd) { return alignmentEnd; }

                    @Override
                    byte[] fetchRef(final RefSequenceSource src, final String chr, final int boundary, final int budget)
                    {
                        return src.getBases(chr, boundary + 1, boundary + budget);
                    }

                    @Override
                    int walk(final byte[] readBases, final byte[] refBases, final int softclipLength, final int walkLength)
                    {
                        // trailing softclip's M boundary is at its start, already boundary-outward
                        final int readStart = readBases.length - softclipLength;
                        final byte[] readWin = Arrays.copyOfRange(readBases, readStart, readStart + walkLength);
                        final byte[] refWin = Arrays.copyOfRange(refBases, 0, walkLength);
                        return BoundaryReclaim.maxScoringPrefix(readWin, refWin);
                    }

                    @Override
                    List<ChrBaseRegion> junctionLookup(
                            final AnnotatedJunctionIndex idx, final String chr, final int boundary, final int offset)
                    {
                        return idx.introByStart(chr, boundary + offset);
                    }
                };

        abstract int softclipIndex(List<CigarElement> cigar);
        abstract int matchedIndex(List<CigarElement> cigar);
        abstract int refBoundary(int alignmentStart, int alignmentEnd);
        abstract byte[] fetchRef(RefSequenceSource src, String chr, int boundary, int budget);
        abstract int walk(byte[] readBases, byte[] refBases, int softclipLength, int walkLength);
        abstract List<ChrBaseRegion> junctionLookup(AnnotatedJunctionIndex idx, String chr, int boundary, int offset);
    }
}
