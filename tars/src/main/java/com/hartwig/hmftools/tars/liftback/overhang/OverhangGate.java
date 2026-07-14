package com.hartwig.hmftools.tars.liftback.overhang;

import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_LENGTH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MIN_OVERHANG_SCORE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Gates a lifted alignment's terminal splice-junction anchors ("overhangs") against the reference genome, and
// reclaims genuine ref-matching bases bwa soft-clipped at read tails.
//
//  peel() - iterative terminal collapse. For each end (trailing then leading) it walks the terminal shape
//     "<...>M nN (anchor)M [(softclip)S]" (or its mirror) inward. It engages a terminal junction only when the trigger applies:
//     an adjacent soft-clip (short overhang with a clip) OR the alignment carries more than one junction (bwa's
//     multi-splice route may not weigh the intronic alternative). A single clean junction with no clip is left
//     untouched. The anchor is scored at its far-exon placement; the junction is kept iff that score exceeds
//     MIN_OVERHANG_SCORE (aligning the anchor beats clipping it), otherwise the N is dropped, the anchor+softclip window is
//     read through the intron contiguously reclaiming its best-scoring prefix to M, and the exposed inner junction
//     becomes terminal and is re-gated. peel() reports dropped=true when it collapses a junction-bearing alignment
//     to a purely contiguous one (no N remains): the caller drops such an XA alt from the tag, but never self.
//
//  reclaimTerminalSoftClip() - the standalone reclaim. For any surviving terminal soft-clip it walks the clip into
//     contiguous reference and converts the best-scoring prefix S -> M. The caller runs this post-supplementary-resolve on
//     the chosen primary only when it is a tx-match, so a genomic/ref over-clip is left as bwa placed it and a clip
//     supplementary resolve needs to splice is never consumed before it runs.
public class OverhangGate
{
    private final RefSequenceSource mRefSource;
    private final OverhangGateStatistics mStatistics;

    public OverhangGate(final RefSequenceSource refSource)
    {
        mRefSource = refSource;
        mStatistics = new OverhangGateStatistics();
    }

    public OverhangGateStatistics statistics()
    {
        return mStatistics;
    }

    // The gate is a no-op without a ref-genome handle (each pass self-guards on it).
    public boolean enabled()
    {
        return mRefSource != null;
    }

    // Final (start, cigar) of a gate pass, plus whether an XA-alt placement should be dropped from the tag.
    public record Result(int pos, String cigar, boolean dropped)
    {
    }

    public Result peel(final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        int pos = alignmentStart;
        String working = cigar;

        // Trailing then leading; each end is collapsed from the outside in, each collapse exposing the next inner junction.
        while(working.indexOf('N') >= 0)
        {
            CollapseResult collapse = tryCollapseTrailing(chromosome, pos, working, readBases);
            if(!collapse.collapsed())
            {
                break;
            }
            working = collapse.cigar();
            mStatistics.countCollapsedTrailing();
        }

        while(working.indexOf('N') >= 0)
        {
            CollapseResult collapse = tryCollapseLeading(chromosome, pos, working, readBases);
            if(!collapse.collapsed())
            {
                break;
            }
            pos = collapse.pos();
            working = collapse.cigar();
            mStatistics.countCollapsedLeading();
        }

        // A junction-bearing alignment collapsed to a purely contiguous one has lost the splice that distinguished
        // it; the caller drops such an XA alt (never self).
        boolean dropped = cigar.indexOf('N') >= 0 && working.indexOf('N') < 0;
        return new Result(pos, working, dropped);
    }

    private CollapseResult tryCollapseTrailing(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
        {
            return CollapseResult.unchanged();
        }

        int last = elements.size() - 1;

        // shape: "<...>M nN yM [eS]"
        boolean hasSoftclip = elements.get(last).getOperator() == CigarOperator.S;
        int softclipLength = hasSoftclip ? elements.get(last).getLength() : 0;
        int anchorIndex = hasSoftclip ? last - 1 : last;
        int intronIndex = anchorIndex - 1;
        int nearIndex = intronIndex - 1;
        if(nearIndex < 0)
        {
            return CollapseResult.unchanged();
        }

        CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()))
        {
            return CollapseResult.unchanged();
        }
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return CollapseResult.unchanged();
        }
        CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
        {
            return CollapseResult.unchanged();
        }

        // trigger: engage only with an adjacent soft-clip or more than one junction (Case 1 / Case 2).
        if(softclipLength == 0 && countIntrons(elements) <= 1)
        {
            return CollapseResult.unchanged();
        }

        int anchorLength = anchor.getLength();
        // a long overhang is trusted outright (a real junction even with a couple of mismatches); only short overhangs
        // are scored. See TarsConstants.MIN_OVERHANG_LENGTH.
        if(anchorLength > MIN_OVERHANG_LENGTH)
        {
            return CollapseResult.unchanged();
        }
        int window = anchorLength + softclipLength;
        int intronLen = elements.get(intronIndex).getLength();
        // genomic end of the near exon, just before the intron
        int nearEnd = alignmentStart + CigarUtils.cigarAlignedLength(elements.subList(0, intronIndex)) - 1;
        byte[] ref = mRefSource.getBases(chromosome, nearEnd + 1, nearEnd + window);
        if(ref == null || ref.length != window)
        {
            return CollapseResult.unchanged();
        }

        // Score the overhang at its far exon (the junction placement). A null far ref (off contig) can't be scored,
        // so it never keeps the junction. Case 1 (soft clip) vs Case 2 (multiple junctions) differ in keepJunction.
        byte[] anchorRead = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length - window + anchorLength);
        byte[] farRef = mRefSource.getBases(chromosome, nearEnd + intronLen + 1, nearEnd + intronLen + anchorLength);
        int overhangScore = (farRef != null && farRef.length == anchorLength)
                ? BoundaryReclaim.score(anchorRead, farRef) : Integer.MIN_VALUE;
        if(keepJunction(overhangScore, anchorRead, ref, softclipLength))
        {
            return CollapseResult.unchanged();
        }

        // drop the intron; reclaim the contiguous-scoring prefix into the near exon, clip the remainder.
        byte[] readWindow = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length);
        int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, ref);
        List<CigarElement> merged = new ArrayList<>(elements.subList(0, nearIndex));
        merged.add(new CigarElement(nearExon.getLength() + reclaimed, CigarOperator.M));
        int residual = window - reclaimed;
        if(residual > 0)
        {
            merged.add(new CigarElement(residual, CigarOperator.S));
        }
        return new CollapseResult(true, alignmentStart, CigarUtils.cigarElementsToStr(merged));
    }

    private CollapseResult tryCollapseLeading(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.size() < 3 || CigarUtils.hasHardClip(elements))
        {
            return CollapseResult.unchanged();
        }

        // shape: "[eS] yM nN <...>M"
        boolean hasSoftclip = elements.get(0).getOperator() == CigarOperator.S;
        int softclipLength = hasSoftclip ? elements.get(0).getLength() : 0;
        int anchorIndex = hasSoftclip ? 1 : 0;
        int intronIndex = anchorIndex + 1;
        int nearIndex = intronIndex + 1;
        if(nearIndex >= elements.size())
        {
            return CollapseResult.unchanged();
        }

        CigarElement anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.getOperator()))
        {
            return CollapseResult.unchanged();
        }
        if(elements.get(intronIndex).getOperator() != CigarOperator.N)
        {
            return CollapseResult.unchanged();
        }
        CigarElement nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.getOperator()))
        {
            return CollapseResult.unchanged();
        }

        // trigger: engage only with an adjacent soft-clip or more than one junction (Case 1 / Case 2).
        if(softclipLength == 0 && countIntrons(elements) <= 1)
        {
            return CollapseResult.unchanged();
        }

        int anchorLength = anchor.getLength();
        // a long overhang is trusted outright (a real junction even with a couple of mismatches); only short overhangs
        // are scored. See TarsConstants.MIN_OVERHANG_LENGTH.
        if(anchorLength > MIN_OVERHANG_LENGTH)
        {
            return CollapseResult.unchanged();
        }
        int window = anchorLength + softclipLength;
        int intron = elements.get(intronIndex).getLength();
        // leading softclip consumes no reference, so nearStart skips only anchor + intron.
        int nearStart = alignmentStart + anchorLength + intron;
        byte[] ref = mRefSource.getBases(chromosome, nearStart - window, nearStart - 1);
        if(ref == null || ref.length != window)
        {
            return CollapseResult.unchanged();
        }

        // the near-exon boundary sits at the END of the leading window, so reverse both so the score walk runs
        // boundary-outward.
        byte[] readWindow = BoundaryReclaim.reversed(Arrays.copyOfRange(readBases, 0, window));
        byte[] nearRef = BoundaryReclaim.reversed(ref);

        // Score the overhang at its far exon (the junction placement); the leading anchor sits at [alignmentStart, +anchorLength).
        // Score is order-independent, so no reversal needed. Case 1 (soft clip) vs Case 2 (multiple junctions) differ
        // in keepJunction.
        byte[] anchorRead = Arrays.copyOfRange(readBases, softclipLength, softclipLength + anchorLength);
        byte[] farRef = mRefSource.getBases(chromosome, alignmentStart, alignmentStart + anchorLength - 1);
        int overhangScore = (farRef != null && farRef.length == anchorLength)
                ? BoundaryReclaim.score(anchorRead, farRef) : Integer.MIN_VALUE;
        if(keepJunction(overhangScore, anchorRead, ref, softclipLength))
        {
            return CollapseResult.unchanged();
        }

        int reclaimed = BoundaryReclaim.maxScoringPrefix(readWindow, nearRef);
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
        return new CollapseResult(true, nearStart - reclaimed, CigarUtils.cigarElementsToStr(merged));
    }

    // Decide whether to keep a short overhang's junction, by case. Case 1 (adjacent soft clip, softclipLength > 0): keep
    // the junction if the overhang's far-exon score beats the clip penalty (MIN_OVERHANG_SCORE), else collapse. Case 2 (more
    // than one junction, no soft clip, softclipLength == 0): keep the junction if the overhang aligns positively (score > 0);
    // only when it does not (score <= 0) is the intronic reference checked - the junction is collapsed only if reading the overhang
    // contiguously on the reference scores strictly higher than the junction (far exon) placement, otherwise it is kept.
    private static boolean keepJunction(final int overhangScore, final byte[] anchorRead, final byte[] windowRef, final int softclipLength)
    {
        if(softclipLength > 0)
        {
            return overhangScore > MIN_OVERHANG_SCORE;
        }
        if(overhangScore > 0)
        {
            return true;
        }
        byte[] contiguousRef = Arrays.copyOfRange(windowRef, 0, anchorRead.length);
        int altRefScore = BoundaryReclaim.score(anchorRead, contiguousRef);
        return altRefScore <= overhangScore;
    }

    // Result of one collapse step. newStart shifts only on a leading collapse.
    private record CollapseResult(boolean collapsed, int pos, String cigar)
    {
        static CollapseResult unchanged()
        {
            return new CollapseResult(false, 0, null);
        }
    }

    public Result reclaimTerminalSoftClip(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.isEmpty() || CigarUtils.hasHardClip(elements))
        {
            return new Result(alignmentStart, cigar, false);
        }

        List<CigarElement> working = new ArrayList<>(elements);
        int alignmentEnd = alignmentStart + CigarUtils.cigarAlignedLength(working) - 1;

        int trailReclaimed = reclaimSide(Side.TRAILING, chromosome, alignmentStart, alignmentEnd, readBases, working);
        int leadReclaimed = reclaimSide(Side.LEADING, chromosome, alignmentStart, alignmentEnd, readBases, working);

        if(trailReclaimed == 0 && leadReclaimed == 0)
        {
            return new Result(alignmentStart, cigar, false);
        }

        mStatistics.countReclaimed(leadReclaimed, trailReclaimed);
        return new Result(alignmentStart - leadReclaimed, CigarUtils.cigarElementsToStr(working), false);
    }

    private int reclaimSide(
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
        if(softclip.getOperator() != CigarOperator.S)
        {
            return 0;
        }
        if(!isMatchedOp(working.get(matchedIndex).getOperator()))
        {
            return 0;
        }

        int boundary = side.refBoundary(alignmentStart, alignmentEnd);

        int budget = softclip.getLength();
        if(side == Side.LEADING && alignmentStart - budget < 1)
        {
            budget = alignmentStart - 1;
        }
        if(budget <= 0)
        {
            return 0;
        }

        byte[] refBases = side.fetchRef(mRefSource, chromosome, boundary, budget);
        if(refBases == null || refBases.length == 0)
        {
            return 0;
        }

        int walkLength = Math.min(refBases.length, budget);
        int reclaimed = side.walk(readBases, refBases, softclip.getLength(), walkLength);
        if(reclaimed == 0)
        {
            return 0;
        }

        applyReclaim(working, softclipIndex, matchedIndex, reclaimed);
        return reclaimed;
    }

    private static void applyReclaim(
            final List<CigarElement> working, final int softclipIndex, final int matchedIndex, final int reclaimed)
    {
        CigarElement softclip = working.get(softclipIndex);
        CigarElement matched = working.get(matchedIndex);

        working.set(matchedIndex, new CigarElement(matched.getLength() + reclaimed, matched.getOperator()));

        if(softclip.getLength() - reclaimed == 0)
        {
            working.remove(softclipIndex);
        }
        else
        {
            working.set(softclipIndex, new CigarElement(softclip.getLength() - reclaimed, CigarOperator.S));
        }
    }

    private static boolean isMatch(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ;
    }

    private static boolean isMatchedOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    private static int countIntrons(final List<CigarElement> elements)
    {
        int count = 0;
        for(CigarElement element : elements)
        {
            if(element.getOperator() == CigarOperator.N)
            {
                ++count;
            }
        }
        return count;
    }

    // Abstracts leading/trailing geometry so reclaimSide is direction-free.
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
                        byte[] readWindow = BoundaryReclaim.reversed(
                                Arrays.copyOfRange(readBases, softclipLength - walkLength, softclipLength));
                        byte[] refWindow = BoundaryReclaim.reversed(
                                Arrays.copyOfRange(refBases, refBases.length - walkLength, refBases.length));
                        return BoundaryReclaim.maxScoringPrefix(readWindow, refWindow);
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
                        byte[] readWindow = Arrays.copyOfRange(readBases, readStart, readStart + walkLength);
                        byte[] refWindow = Arrays.copyOfRange(refBases, 0, walkLength);
                        return BoundaryReclaim.maxScoringPrefix(readWindow, refWindow);
                    }
                };

        abstract int softclipIndex(List<CigarElement> cigar);

        abstract int matchedIndex(List<CigarElement> cigar);

        abstract int refBoundary(int alignmentStart, int alignmentEnd);

        abstract byte[] fetchRef(RefSequenceSource src, String chr, int boundary, int budget);

        abstract int walk(byte[] readBases, byte[] refBases, int softclipLength, int walkLength);
    }
}
