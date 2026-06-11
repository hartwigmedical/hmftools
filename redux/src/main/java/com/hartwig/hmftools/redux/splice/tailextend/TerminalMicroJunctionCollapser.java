package com.hartwig.hmftools.redux.splice.tailextend;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SKIPPED;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SOFTCLIP;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.redux.splice.rescue.CigarShape;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;

// Removes a spurious terminal junction fabricated when the tx-contig walk over-extends a read past an
// exon boundary. The tiny anchor (yM, y <= MaxAnchor) is collapsed when its bases match the genome
// contiguously past the near exon; the intron is dropped and any reclaimed softclip bases are merged
// into the near exon. Anchor bases must match exactly (no error budget); softclip bases tolerate
// max(1, length/10) mismatches. Legitimate short-anchor junctions are preserved because their anchor
// bases do NOT match the contiguous genome.
public class TerminalMicroJunctionCollapser
{
    private final RefSequenceSource mRefSource;
    private final int mMaxAnchor;
    private final AtomicLong mCollapsedLeading = new AtomicLong();
    private final AtomicLong mCollapsedTrailing = new AtomicLong();

    public TerminalMicroJunctionCollapser(final RefSequenceSource refSource, final int maxAnchor)
    {
        mRefSource = refSource;
        mMaxAnchor = maxAnchor;
    }

    public long collapsedLeading() { return mCollapsedLeading.get(); }
    public long collapsedTrailing() { return mCollapsedTrailing.get(); }

    public TerminalCollapseResult tryCollapse(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
            return TerminalCollapseResult.unchanged();

        final List<CigarShape.Element> elements = CigarShape.parse(cigar);
        if(elements.size() < 3 || CigarShape.hasHardClip(elements))
            return TerminalCollapseResult.unchanged();

        // both ends are independent; apply trailing first, then leading on the result.
        final TerminalCollapseResult trailing = tryTrailing(chromosome, alignmentStart, elements, readBases);
        final List<CigarShape.Element> afterTrailing = trailing.Collapsed ? CigarShape.parse(trailing.NewCigar) : elements;
        final int startAfterTrailing = trailing.Collapsed ? trailing.NewStart : alignmentStart;

        final TerminalCollapseResult leading = tryLeading(chromosome, startAfterTrailing, afterTrailing, readBases);
        if(leading.Collapsed)
            return leading;
        return trailing;
    }

    private TerminalCollapseResult tryTrailing(
            final String chromosome, final int alignmentStart, final List<CigarShape.Element> elements, final byte[] readBases)
    {
        final int last = elements.size() - 1;

        // shape: "<...>M nN yM [eS]"
        final boolean hasSoftclip = elements.get(last).Op == OP_SOFTCLIP;
        final int e = hasSoftclip ? elements.get(last).Length : 0;
        final int anchorIndex = hasSoftclip ? last - 1 : last;
        final int intronIndex = anchorIndex - 1;
        final int nearIndex = intronIndex - 1;
        if(nearIndex < 0)
            return TerminalCollapseResult.unchanged();

        final CigarShape.Element anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.Op) || anchor.Length > mMaxAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(intronIndex).Op != OP_SKIPPED)
            return TerminalCollapseResult.unchanged();
        final CigarShape.Element nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.Op))
            return TerminalCollapseResult.unchanged();

        final int y = anchor.Length;
        final int window = y + e;
        // genomic end of the near exon, just before the intron
        final int nearEnd = alignmentStart + CigarShape.referenceSpan(elements.subList(0, intronIndex)) - 1;
        final byte[] ref = mRefSource.getBases(chromosome, nearEnd + 1, nearEnd + window);
        if(ref == null || ref.length != window)
            return TerminalCollapseResult.unchanged();

        // read window measured from the near-exon boundary outward: anchor then softclip
        final byte[] readWindow = Arrays.copyOfRange(readBases, readBases.length - window, readBases.length);
        final int reclaimed = contiguousRun(readWindow, ref, y, e, true);
        if(reclaimed == 0)
            return TerminalCollapseResult.unchanged();

        // drop intron + anchor, grow near exon; keep unmatched softclip remainder.
        final List<CigarShape.Element> merged = new ArrayList<>(elements.subList(0, nearIndex));
        merged.add(new CigarShape.Element(nearExon.Length + reclaimed, OP_MATCH));
        final int residual = window - reclaimed;
        if(residual > 0)
            merged.add(new CigarShape.Element(residual, OP_SOFTCLIP));
        mCollapsedTrailing.incrementAndGet();
        return new TerminalCollapseResult(true, alignmentStart, CigarShape.format(merged));
    }

    private TerminalCollapseResult tryLeading(
            final String chromosome, final int alignmentStart, final List<CigarShape.Element> elements, final byte[] readBases)
    {
        // shape: "[eS] yM nN <...>M"
        final boolean hasSoftclip = elements.get(0).Op == OP_SOFTCLIP;
        final int e = hasSoftclip ? elements.get(0).Length : 0;
        final int anchorIndex = hasSoftclip ? 1 : 0;
        final int intronIndex = anchorIndex + 1;
        final int nearIndex = intronIndex + 1;
        if(nearIndex >= elements.size())
            return TerminalCollapseResult.unchanged();

        final CigarShape.Element anchor = elements.get(anchorIndex);
        if(!isMatch(anchor.Op) || anchor.Length > mMaxAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(intronIndex).Op != OP_SKIPPED)
            return TerminalCollapseResult.unchanged();
        final CigarShape.Element nearExon = elements.get(nearIndex);
        if(!isMatch(nearExon.Op))
            return TerminalCollapseResult.unchanged();

        final int y = anchor.Length;
        final int window = y + e;
        final int intron = elements.get(intronIndex).Length;
        // leading softclip consumes no reference, so nearStart skips only anchor + intron.
        final int nearStart = alignmentStart + y + intron;
        final byte[] ref = mRefSource.getBases(chromosome, nearStart - window, nearStart - 1);
        if(ref == null || ref.length != window)
            return TerminalCollapseResult.unchanged();

        // walk from the near-exon end outward (anchor first, then softclip).
        final byte[] readWindow = Arrays.copyOfRange(readBases, 0, window);
        final int reclaimed = contiguousRun(readWindow, ref, y, e, false);
        if(reclaimed == 0)
            return TerminalCollapseResult.unchanged();

        final int residual = window - reclaimed;
        final List<CigarShape.Element> merged = new ArrayList<>();
        if(residual > 0)
            merged.add(new CigarShape.Element(residual, OP_SOFTCLIP));
        merged.add(new CigarShape.Element(reclaimed + nearExon.Length, OP_MATCH));
        for(int i = nearIndex + 1; i < elements.size(); ++i)
            merged.add(elements.get(i));
        mCollapsedLeading.incrementAndGet();
        return new TerminalCollapseResult(true, nearStart - reclaimed, CigarShape.format(merged));
    }

    // Returns reclaimed length in [y, y+e]: anchor (y bases, must match exactly), then softclip
    // (tolerates max(1, length/10) mismatches). anchorAtStart=true for trailing, false for leading.
    private static int contiguousRun(
            final byte[] readWindow, final byte[] ref, final int y, final int e, final boolean anchorAtStart)
    {
        final int window = y + e;
        int mismatches = 0;
        int reclaimed = 0;
        for(int step = 0; step < window; ++step)
        {
            final int idx = anchorAtStart ? step : window - 1 - step;
            final boolean match = basesEqualIgnoreCase(readWindow[idx], ref[idx]);
            final int length = step + 1;
            if(step < y)
            {
                if(!match)
                    return 0;
                reclaimed = length;
                continue;
            }
            if(!match)
                ++mismatches;
            final int allowed = Math.max(1, length / 10);
            if(mismatches <= allowed)
                reclaimed = length;
            else if(mismatches > allowed + 1)
                break;
        }
        return reclaimed;
    }

    private static boolean isMatch(final char op)
    {
        return op == OP_MATCH || op == OP_SEQ_MATCH;
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        return Character.toUpperCase((char) a) == Character.toUpperCase((char) b);
    }
}
