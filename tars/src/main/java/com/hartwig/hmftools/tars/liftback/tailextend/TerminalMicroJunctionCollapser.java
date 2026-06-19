package com.hartwig.hmftools.tars.liftback.tailextend;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.tars.common.SpliceCommon;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Resolves an untrusted terminal junction fabricated when the tx-contig walk over-extends a read past an
// exon boundary, shape "<...>M nN yM [eS]" (trailing) or "[eS] yM nN <...>M" (leading) where the terminal
// anchor y is below SpliceCommon.MIN_JUNCTION_ANCHOR. A sub-threshold anchor carries too little evidence to
// assert a junction on its own, so liftback re-evaluates it rather than trusting bwa's placement.
//
// Decision is a head-to-head alignment-score comparison: score bwa's split (the anchor on the far exon, across
// the intron; clip bases discarded) against the whole terminal window extended contiguously into the intron
// (no junction). If the split scores strictly higher it is a real short-overhang junction and is kept;
// otherwise the contiguous extension explains the tail at least as well, the intron is dropped, the
// contiguous-matching prefix is reclaimed into the near exon as M and the rest is soft-clipped (e.g.
// "...M xN 2M 1S" -> "...M 3S"; intron retention reclaims the whole tail).
public class TerminalMicroJunctionCollapser
{
    private final RefSequenceSource mRefSource;
    private final int mMinJunctionAnchor;
    private long mCollapsedLeading;
    private long mCollapsedTrailing;

    public TerminalMicroJunctionCollapser(final RefSequenceSource refSource, final int minJunctionAnchor)
    {
        mRefSource = refSource;
        mMinJunctionAnchor = minJunctionAnchor;
    }

    public long collapsedLeading() { return mCollapsedLeading; }

    public long collapsedTrailing() { return mCollapsedTrailing; }

    public TerminalCollapseResult tryCollapse(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
            return TerminalCollapseResult.unchanged();

        final List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
            return TerminalCollapseResult.unchanged();

        // both ends are independent; apply trailing first, then leading on the result.
        final TerminalCollapseResult trailing = tryTrailing(chromosome, alignmentStart, elements, readBases);
        final List<CigarElement> afterTrailing = trailing.collapsed() ? CigarUtils.cigarElementsFromStr(trailing.newCigar()) : elements;
        final int startAfterTrailing = trailing.collapsed() ? trailing.newStart() : alignmentStart;

        final TerminalCollapseResult leading = tryLeading(chromosome, startAfterTrailing, afterTrailing, readBases);
        if(leading.collapsed())
            return leading;
        return trailing;
    }

    private TerminalCollapseResult tryTrailing(
            final String chromosome, final int alignmentStart, final List<CigarElement> elements, final byte[] readBases)
    {
        final int last = elements.size() - 1;

        // shape: "<...>M nN yM [eS]"
        final boolean hasSoftclip = elements.get(last).getOperator() == CigarOperator.S;
        final int e = hasSoftclip ? elements.get(last).getLength() : 0;
        final int anchorIndex = hasSoftclip ? last - 1 : last;
        final int intronIndex = anchorIndex - 1;
        final int nearIndex = intronIndex - 1;
        if(nearIndex < 0)
            return TerminalCollapseResult.unchanged();

        final CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()) || anchor.getLength() >= mMinJunctionAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
            return TerminalCollapseResult.unchanged();
        final CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
            return TerminalCollapseResult.unchanged();

        final int y = anchor.getLength();
        final int window = y + e;
        final int intronLen = elements.get(intronIndex).getLength();
        // genomic end of the near exon, just before the intron
        final int nearEnd = alignmentStart + CigarUtils.cigarAlignedLength(elements.subList(0, intronIndex)) - 1;
        final byte[] ref = mRefSource.getBases(chromosome, nearEnd + 1, nearEnd + window);
        if(ref == null || ref.length != window)
            return TerminalCollapseResult.unchanged();

        // re-score the terminal region against the contiguous near-exon genome; the boundary sits at the
        // start of the window, so walk forward.
        final byte[] readWindow = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length);
        final int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, ref);

        // Keep bwa's junction only when the anchor scores strictly better on the far exon (bwa's placement,
        // across the intron) than the whole terminal window does extended contiguously into the intron. That is
        // a real short-overhang junction. Otherwise the contiguous extension explains the tail at least as well,
        // so the sub-threshold intron is a tx over-run and is dropped. Clip bases score 0 in the split (clipped),
        // so a tail that only matches contiguously (intron retention) always loses to extension here.
        final byte[] anchorRead = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length - window + y);
        final byte[] farRef = mRefSource.getBases(chromosome, nearEnd + intronLen + 1, nearEnd + intronLen + y);
        if(splitWins(anchorRead, farRef, BoundaryReclaim.maxPrefixScore(readWindow, ref)))
            return TerminalCollapseResult.unchanged();

        // drop the intron; reclaim the scoring prefix into the near exon, clip the remainder.
        final List<CigarElement> merged = new ArrayList<>(elements.subList(0, nearIndex));
        merged.add(new CigarElement(nearExon.getLength() + reclaimed, CigarOperator.M));
        final int residual = window - reclaimed;
        if(residual > 0)
            merged.add(new CigarElement(residual, CigarOperator.S));
        ++mCollapsedTrailing;
        return new TerminalCollapseResult(true, alignmentStart, CigarUtils.cigarElementsToStr(merged));
    }

    private TerminalCollapseResult tryLeading(
            final String chromosome, final int alignmentStart, final List<CigarElement> elements, final byte[] readBases)
    {
        // shape: "[eS] yM nN <...>M"
        final boolean hasSoftclip = elements.get(0).getOperator() == CigarOperator.S;
        final int e = hasSoftclip ? elements.get(0).getLength() : 0;
        final int anchorIndex = hasSoftclip ? 1 : 0;
        final int intronIndex = anchorIndex + 1;
        final int nearIndex = intronIndex + 1;
        if(nearIndex >= elements.size())
            return TerminalCollapseResult.unchanged();

        final CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()) || anchor.getLength() >= mMinJunctionAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
            return TerminalCollapseResult.unchanged();
        final CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
            return TerminalCollapseResult.unchanged();

        final int y = anchor.getLength();
        final int window = y + e;
        final int intron = elements.get(intronIndex).getLength();
        // leading softclip consumes no reference, so nearStart skips only anchor + intron.
        final int nearStart = alignmentStart + y + intron;
        final byte[] ref = mRefSource.getBases(chromosome, nearStart - window, nearStart - 1);
        if(ref == null || ref.length != window)
            return TerminalCollapseResult.unchanged();

        // the near-exon boundary sits at the END of the leading window, so reverse both so the score walk
        // runs boundary-outward.
        final byte[] readWindow = BoundaryReclaim.reversed(Arrays.copyOfRange(readBases, 0, window));
        final byte[] nearRef = BoundaryReclaim.reversed(ref);
        final int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, nearRef);

        // keep bwa's junction only when the anchor scores strictly better on the far exon than the contiguous
        // extension does (see tryTrailing). The leading anchor sits at [alignmentStart, +y) - the read's leading
        // bases after the softclip. Split score is order-independent, so no reversal needed.
        final byte[] anchorRead = Arrays.copyOfRange(readBases, e, e + y);
        final byte[] farRef = mRefSource.getBases(chromosome, alignmentStart, alignmentStart + y - 1);
        if(splitWins(anchorRead, farRef, BoundaryReclaim.maxPrefixScore(readWindow, nearRef)))
            return TerminalCollapseResult.unchanged();

        final int residual = window - reclaimed;
        final List<CigarElement> merged = new ArrayList<>();
        if(residual > 0)
            merged.add(new CigarElement(residual, CigarOperator.S));
        merged.add(new CigarElement(reclaimed + nearExon.getLength(), CigarOperator.M));
        for(int i = nearIndex + 1; i < elements.size(); ++i)
            merged.add(elements.get(i));
        ++mCollapsedLeading;
        return new TerminalCollapseResult(true, nearStart - reclaimed, CigarUtils.cigarElementsToStr(merged));
    }

    // bwa's split placement wins when the anchor scores strictly better on the far exon than the whole terminal
    // window scores extended contiguously into the intron. A null far ref (off contig) can't beat anything, so
    // it never blocks the collapse.
    private static boolean splitWins(final byte[] anchorRead, final byte[] farRef, final int contiguousScore)
    {
        if(farRef == null || farRef.length != anchorRead.length)
            return false;
        return BoundaryReclaim.score(anchorRead, farRef) > contiguousScore;
    }

    private static boolean isMatch(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ;
    }
}
