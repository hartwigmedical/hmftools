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

// Reconciles a lifted placement against the genome with two terminal passes that share a ref-genome handle:
//
//  1. terminal micro-junction collapse - resolves an untrusted terminal junction fabricated when the tx-contig
//     walk over-extends a read past an exon boundary, shape "<...>M nN yM [eS]" (trailing) or "[eS] yM nN <...>M"
//     (leading) where the terminal anchor y is below TarsConstants.MIN_JUNCTION_ANCHOR. A sub-threshold anchor
//     carries too little evidence to assert a junction on its own, so liftback re-evaluates it rather than
//     trusting bwa's placement. Decision is a head-to-head alignment-score comparison: score bwa's split (the
//     anchor on the far exon, across the intron; clip bases discarded) against the whole terminal window extended
//     contiguously into the intron (no junction). If the split scores strictly higher it is a real short-overhang
//     junction and is kept; otherwise the contiguous extension explains the tail at least as well, the intron is
//     dropped, the contiguous-matching prefix is reclaimed into the near exon as M and the rest is soft-clipped
//     (e.g. "...M xN 2M 1S" -> "...M 3S"; intron retention reclaims the whole tail).
//
//  2. softclip tail extension - recovers ref-matching bases bwa-mem2 softclipped at read tails. Runs after
//     junction rescue so rescue's boundary lookups see the original cigar. Skips when an annotated intron starts
//     within MaxExtension of the softclip boundary - that range belongs to junction rescue.
public class TerminalReconciler
{
    private final RefSequenceSource mRefSource;
    private final int mMinJunctionAnchor;
    private final AnnotatedJunctionIndex mJunctionGuard;
    private final TailExtensionConfig mConfig;
    private final TailExtensionStatistics mStatistics;

    private long mCollapsedLeading;
    private long mCollapsedTrailing;

    public TerminalReconciler(
            final RefSequenceSource refSource, final int minJunctionAnchor,
            final AnnotatedJunctionIndex junctionGuard, final TailExtensionConfig config)
    {
        mRefSource = refSource;
        mMinJunctionAnchor = minJunctionAnchor;
        mJunctionGuard = junctionGuard;
        mConfig = config;
        mStatistics = new TailExtensionStatistics();
    }

    public TailExtensionStatistics statistics()
    {
        return mStatistics;
    }

    public long collapsedLeading() { return mCollapsedLeading; }

    public long collapsedTrailing() { return mCollapsedTrailing; }

    // Restore the collapse counters to roll back a discarded provisional mate decision without double-counting
    // (see LiftBackGroupProcessor.processNameGroup).
    public void restoreCollapseCounters(final long collapsedLeading, final long collapsedTrailing)
    {
        mCollapsedLeading = collapsedLeading;
        mCollapsedTrailing = collapsedTrailing;
    }

    // Collapse is a no-op without a ref-genome handle (tryCollapse self-guards on it). Mirrors the old
    // "terminal collapser constructed only when refSource != null" gating in the orchestrator.
    public boolean collapseEnabled() { return mRefSource != null; }

    // Reconcile one lifted placement: collapse a fabricated terminal micro-junction first (only when the cigar has
    // an N), then walk a reclaimable terminal softclip into contiguous genome (only when allowTailExtend and the
    // cigar has an S). Replicates the exact pass sequence the orchestrator used to run inline per candidate.
    public ReconcileResult reconcile(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases,
            final boolean allowTailExtend)
    {
        int pos = alignmentStart;
        String workingCigar = cigar;

        if(workingCigar != null && workingCigar.indexOf('N') >= 0)
        {
            TerminalCollapseResult collapse = tryCollapse(chromosome, pos, workingCigar, readBases);
            if(collapse.collapsed())
            {
                pos = collapse.newStart();
                workingCigar = collapse.newCigar();
            }
        }

        if(allowTailExtend && workingCigar != null && workingCigar.indexOf('S') >= 0)
        {
            TailExtensionResult extension = tryExtend(chromosome, pos, workingCigar, readBases);
            if(extension.extended())
            {
                pos = extension.newStart();
                workingCigar = extension.newCigar();
            }
        }

        return new ReconcileResult(pos, workingCigar);
    }

    // Final (start, cigar) of a reconcile pass.
    public record ReconcileResult(int pos, String cigar)
    {
    }

    // ============================ terminal micro-junction collapse ============================

    public TerminalCollapseResult tryCollapse(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return TerminalCollapseResult.unchanged();
        }

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
        {
            return TerminalCollapseResult.unchanged();
        }

        // both ends are independent; apply trailing first, then leading on the result.
        TerminalCollapseResult trailing = tryTrailing(chromosome, alignmentStart, elements, readBases);
        List<CigarElement> afterTrailing = trailing.collapsed() ? CigarUtils.cigarElementsFromStr(trailing.newCigar()) : elements;
        int startAfterTrailing = trailing.collapsed() ? trailing.newStart() : alignmentStart;

        TerminalCollapseResult leading = tryLeading(chromosome, startAfterTrailing, afterTrailing, readBases);
        if(leading.collapsed())
        {
            return leading;
        }
        return trailing;
    }

    private TerminalCollapseResult tryTrailing(
            final String chromosome, final int alignmentStart, final List<CigarElement> elements, final byte[] readBases)
    {
        int last = elements.size() - 1;

        // shape: "<...>M nN yM [eS]"
        boolean hasSoftclip = elements.get(last).getOperator() == CigarOperator.S;
        int e = hasSoftclip ? elements.get(last).getLength() : 0;
        int anchorIndex = hasSoftclip ? last - 1 : last;
        int intronIndex = anchorIndex - 1;
        int nearIndex = intronIndex - 1;
        if(nearIndex < 0)
        {
            return TerminalCollapseResult.unchanged();
        }

        CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()) || anchor.getLength() >= mMinJunctionAnchor)
        {
            return TerminalCollapseResult.unchanged();
        }
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return TerminalCollapseResult.unchanged();
        }
        CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
        {
            return TerminalCollapseResult.unchanged();
        }

        int y = anchor.getLength();
        int window = y + e;
        int intronLen = elements.get(intronIndex).getLength();
        // genomic end of the near exon, just before the intron
        int nearEnd = alignmentStart + CigarUtils.cigarAlignedLength(elements.subList(0, intronIndex)) - 1;
        byte[] ref = mRefSource.getBases(chromosome, nearEnd + 1, nearEnd + window);
        if(ref == null || ref.length != window)
        {
            return TerminalCollapseResult.unchanged();
        }

        // re-score the terminal region against the contiguous near-exon genome; the boundary sits at the
        // start of the window, so walk forward.
        byte[] readWindow = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length);
        int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, ref);

        // Keep bwa's junction only when the anchor scores strictly better on the far exon (bwa's placement,
        // across the intron) than the whole terminal window does extended contiguously into the intron. That is
        // a real short-overhang junction. Otherwise the contiguous extension explains the tail at least as well,
        // so the sub-threshold intron is a tx over-run and is dropped. Clip bases score 0 in the split (clipped),
        // so a tail that only matches contiguously (intron retention) always loses to extension here.
        byte[] anchorRead = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length - window + y);
        byte[] farRef = mRefSource.getBases(chromosome, nearEnd + intronLen + 1, nearEnd + intronLen + y);
        if(splitWins(anchorRead, farRef, BoundaryReclaim.maxPrefixScore(readWindow, ref)))
        {
            return TerminalCollapseResult.unchanged();
        }

        // drop the intron; reclaim the scoring prefix into the near exon, clip the remainder.
        List<CigarElement> merged = new ArrayList<>(elements.subList(0, nearIndex));
        merged.add(new CigarElement(nearExon.getLength() + reclaimed, CigarOperator.M));
        int residual = window - reclaimed;
        if(residual > 0)
        {
            merged.add(new CigarElement(residual, CigarOperator.S));
        }
        ++mCollapsedTrailing;
        return new TerminalCollapseResult(true, alignmentStart, CigarUtils.cigarElementsToStr(merged));
    }

    private TerminalCollapseResult tryLeading(
            final String chromosome, final int alignmentStart, final List<CigarElement> elements, final byte[] readBases)
    {
        // shape: "[eS] yM nN <...>M"
        boolean hasSoftclip = elements.get(0).getOperator() == CigarOperator.S;
        int e = hasSoftclip ? elements.get(0).getLength() : 0;
        int anchorIndex = hasSoftclip ? 1 : 0;
        int intronIndex = anchorIndex + 1;
        int nearIndex = intronIndex + 1;
        if(nearIndex >= elements.size())
        {
            return TerminalCollapseResult.unchanged();
        }

        CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()) || anchor.getLength() >= mMinJunctionAnchor)
        {
            return TerminalCollapseResult.unchanged();
        }
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return TerminalCollapseResult.unchanged();
        }
        CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
        {
            return TerminalCollapseResult.unchanged();
        }

        int y = anchor.getLength();
        int window = y + e;
        int intron = elements.get(intronIndex).getLength();
        // leading softclip consumes no reference, so nearStart skips only anchor + intron.
        int nearStart = alignmentStart + y + intron;
        byte[] ref = mRefSource.getBases(chromosome, nearStart - window, nearStart - 1);
        if(ref == null || ref.length != window)
        {
            return TerminalCollapseResult.unchanged();
        }

        // the near-exon boundary sits at the END of the leading window, so reverse both so the score walk
        // runs boundary-outward.
        byte[] readWindow = BoundaryReclaim.reversed(Arrays.copyOfRange(readBases, 0, window));
        byte[] nearRef = BoundaryReclaim.reversed(ref);
        int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, nearRef);

        // keep bwa's junction only when the anchor scores strictly better on the far exon than the contiguous
        // extension does (see tryTrailing). The leading anchor sits at [alignmentStart, +y) - the read's leading
        // bases after the softclip. Split score is order-independent, so no reversal needed.
        byte[] anchorRead = Arrays.copyOfRange(readBases, e, e + y);
        byte[] farRef = mRefSource.getBases(chromosome, alignmentStart, alignmentStart + y - 1);
        if(splitWins(anchorRead, farRef, BoundaryReclaim.maxPrefixScore(readWindow, nearRef)))
        {
            return TerminalCollapseResult.unchanged();
        }

        int residual = window - reclaimed;
        List<CigarElement> merged = new ArrayList<>();
        if(residual > 0)
        {
            merged.add(new CigarElement(residual, CigarOperator.S));
        }
        merged.add(new CigarElement(reclaimed + nearExon.getLength(), CigarOperator.M));
        for(int i = nearIndex + 1; i < elements.size(); ++i)
        {
            merged.add(elements.get(i));
        }
        ++mCollapsedLeading;
        return new TerminalCollapseResult(true, nearStart - reclaimed, CigarUtils.cigarElementsToStr(merged));
    }

    // bwa's split placement wins only when the anchor scores strictly better on the far exon than the whole
    // terminal window scores extended contiguously into the intron. A real short junction's anchor matches the
    // far exon but not the contiguous near genome, so it wins strictly and is kept; a fabricated anchor matches
    // both (a tie) or the contiguous walk better, so it collapses. A null far ref (off contig) can't beat
    // anything, so it never blocks the collapse.
    private static boolean splitWins(final byte[] anchorRead, final byte[] farRef, final int contiguousScore)
    {
        if(farRef == null || farRef.length != anchorRead.length)
        {
            return false;
        }
        return BoundaryReclaim.score(anchorRead, farRef) > contiguousScore;
    }

    private static boolean isMatch(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ;
    }

    // ============================ softclip tail extension ============================

    public TailExtensionResult tryExtend(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(!mConfig.Enabled || mRefSource == null || chromosome == null
                || cigar == null || readBases == null || readBases.length == 0)
            return TailExtensionResult.unchanged();

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.isEmpty())
        {
            return TailExtensionResult.unchanged();
        }
        if(CigarUtils.hasHardClip(elements))
        {
            mStatistics.countSkippedComplexShape();
            return TailExtensionResult.unchanged();
        }

        mStatistics.countEvaluated();

        List<CigarElement> working = new ArrayList<>(elements);
        int alignmentEnd = alignmentStart + CigarUtils.cigarAlignedLength(working) - 1;

        int trailExtended = tryExtendSide(Side.TRAILING, chromosome, alignmentStart, alignmentEnd, readBases, working);
        int leadExtended = tryExtendSide(Side.LEADING, chromosome, alignmentStart, alignmentEnd, readBases, working);

        if(trailExtended == 0 && leadExtended == 0)
        {
            return TailExtensionResult.unchanged();
        }

        mStatistics.countExtended(leadExtended, trailExtended);
        int newStart = alignmentStart - leadExtended;
        return new TailExtensionResult(true, newStart, CigarUtils.cigarElementsToStr(working),
                leadExtended, trailExtended);
    }

    private int tryExtendSide(
            final Side side, final String chromosome, final int alignmentStart, final int alignmentEnd,
            final byte[] readBases, final List<CigarElement> working)
    {
        if(working.size() < 2)
        {
            return 0;
        }

        int softclipIndex = side.softclipIndex(working);
        int matchedIndex = side.matchedIndex(working);

        CigarElement softclip = working.get(softclipIndex);
        if(softclip.getOperator() != CigarOperator.S || softclip.getLength() < mConfig.MinSoftclipLength)
        {
            return 0;
        }

        if(!isMatchedOp(working.get(matchedIndex).getOperator()))
        {
            mStatistics.countSkippedComplexShape();
            return 0;
        }

        int boundary = side.refBoundary(alignmentStart, alignmentEnd);

        int extendBudget = Math.min(softclip.getLength(), mConfig.MaxExtension);
        if(side == Side.LEADING && alignmentStart - extendBudget < 1)
        {
            extendBudget = alignmentStart - 1;
        }
        if(extendBudget < mConfig.MinExtension)
        {
            return 0;
        }

        byte[] refBases = side.fetchRef(mRefSource, chromosome, boundary, extendBudget);
        if(refBases == null || refBases.length < mConfig.MinExtension)
        {
            mStatistics.countSkippedNoRef();
            return 0;
        }

        int walkLength = Math.min(refBases.length, extendBudget);
        int bestLength = side.walk(readBases, refBases, softclip.getLength(), walkLength);

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
            {
                mStatistics.countRejectedTooManyMismatches();
            }
            return 0;
        }

        applyExtension(working, softclipIndex, matchedIndex, bestLength);
        return bestLength;
    }

    private boolean nearAnnotatedJunction(final Side side, final String chromosome, final int boundary)
    {
        if(mJunctionGuard == null)
        {
            return false;
        }
        for(int offset = 1; offset <= mConfig.MaxExtension; ++offset)
        {
            if(!side.junctionLookup(mJunctionGuard, chromosome, boundary, offset).isEmpty())
            {
                return true;
            }
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
        CigarElement softclip = working.get(softclipIndex);
        CigarElement matched = working.get(matchedIndex);

        CigarElement grownMatch = new CigarElement(matched.getLength() + extend, matched.getOperator());
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
                        byte[] readWin = BoundaryReclaim.reversed(
                                Arrays.copyOfRange(readBases, softclipLength - walkLength, softclipLength));
                        byte[] refWin = BoundaryReclaim.reversed(
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
                        int readStart = readBases.length - softclipLength;
                        byte[] readWin = Arrays.copyOfRange(readBases, readStart, readStart + walkLength);
                        byte[] refWin = Arrays.copyOfRange(refBases, 0, walkLength);
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
