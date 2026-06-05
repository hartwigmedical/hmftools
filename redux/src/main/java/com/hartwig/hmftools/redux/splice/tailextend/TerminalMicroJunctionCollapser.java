package com.hartwig.hmftools.redux.splice.tailextend;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SKIPPED;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.redux.splice.rescue.CigarShape;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;

// Removes a spurious junction that the tx-contig walk fabricated at a read terminus. When a read's
// bases run a base or two past an exon boundary on the concatenated transcript contig, lift-back
// injects a tiny terminal anchor across an intron: "<...>M nN yM" (trailing) or "yM nN <...>M"
// (leading), with y at or below MaxAnchor. STAR routinely aligns these same bases contiguously (no
// intron), because the bases also match the genome continuing the near exon. When they do, the
// junction is unnecessary and is collapsed: the yM is merged into the adjacent exon's M and the nN is
// dropped (a leading collapse moves the alignment start back by y+n).
//
// This is NOT the legitimate short-anchor case (e.g. a real "2M 80N 149M" where the read genuinely
// starts 2 bases into an exon): there the terminal bases do NOT match the genome contiguously (an
// intron sits between), so the ref check fails and the junction is kept. The ref comparison is the
// discriminator, so anchor length alone never drives the decision.
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

        // try trailing first, then leading; a read can only have one terminal micro-junction per end so
        // at most one collapse per side, but both ends are independent.
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
        final CigarShape.Element tail = elements.get(last);
        if(!isMatch(tail.Op) || tail.Length > mMaxAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(last - 1).Op != OP_SKIPPED)
            return TerminalCollapseResult.unchanged();
        final CigarShape.Element nearExon = elements.get(last - 2);
        if(!isMatch(nearExon.Op))
            return TerminalCollapseResult.unchanged();

        final int y = tail.Length;
        // genomic end of the near exon's matched run (before the intron)
        final int nearEnd = alignmentStart + CigarShape.referenceSpan(elements.subList(0, last - 1)) - 1;
        final byte[] ref = mRefSource.getBases(chromosome, nearEnd + 1, nearEnd + y);
        if(ref == null || ref.length != y)
            return TerminalCollapseResult.unchanged();

        // the trailing y read bases (no trailing softclip in this shape) align contiguously?
        if(!tailReadBasesMatchRef(readBases, readBases.length - y, ref))
            return TerminalCollapseResult.unchanged();

        // collapse: drop the intron + tiny anchor, extend the near exon by y
        final List<CigarShape.Element> merged = new ArrayList<>(elements.subList(0, last - 2));
        merged.add(new CigarShape.Element(nearExon.Length + y, OP_MATCH));
        mCollapsedTrailing.incrementAndGet();
        return new TerminalCollapseResult(true, alignmentStart, CigarShape.format(merged));
    }

    private TerminalCollapseResult tryLeading(
            final String chromosome, final int alignmentStart, final List<CigarShape.Element> elements, final byte[] readBases)
    {
        final CigarShape.Element head = elements.get(0);
        if(!isMatch(head.Op) || head.Length > mMaxAnchor)
            return TerminalCollapseResult.unchanged();
        if(elements.get(1).Op != OP_SKIPPED)
            return TerminalCollapseResult.unchanged();
        final CigarShape.Element nearExon = elements.get(2);
        if(!isMatch(nearExon.Op))
            return TerminalCollapseResult.unchanged();

        final int y = head.Length;
        final int intron = elements.get(1).Length;
        // genomic start of the near exon's matched run (after the intron)
        final int nearStart = alignmentStart + y + intron;
        final byte[] ref = mRefSource.getBases(chromosome, nearStart - y, nearStart - 1);
        if(ref == null || ref.length != y)
            return TerminalCollapseResult.unchanged();

        // the leading y read bases align contiguously continuing the near exon backwards?
        if(!tailReadBasesMatchRef(readBases, 0, ref))
            return TerminalCollapseResult.unchanged();

        final List<CigarShape.Element> merged = new ArrayList<>();
        merged.add(new CigarShape.Element(y + nearExon.Length, OP_MATCH));
        for(int i = 3; i < elements.size(); ++i)
            merged.add(elements.get(i));
        mCollapsedLeading.incrementAndGet();
        return new TerminalCollapseResult(true, nearStart - y, CigarShape.format(merged));
    }

    // exact match required: a 1-2bp terminal anchor carries no error budget, and the whole point is
    // that the bases coincide with the contiguous genome.
    private static boolean tailReadBasesMatchRef(final byte[] readBases, final int readOffset, final byte[] ref)
    {
        for(int i = 0; i < ref.length; ++i)
        {
            if(!basesEqualIgnoreCase(readBases[readOffset + i], ref[i]))
                return false;
        }
        return true;
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
