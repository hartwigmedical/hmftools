package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SEQ_MISMATCH;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.OP_SKIPPED;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.consumesRead;
import static com.hartwig.hmftools.redux.splice.rescue.CigarShape.consumesReference;

import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.redux.splice.rescue.CigarShape;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;
import com.hartwig.hmftools.redux.splice.rescue.SpliceMotif;

// Slides a lifted intron a few bp to land its donor/acceptor on a canonical splice motif (GT-AG, or
// the reverse-strand complement CT-AC). A tx-contig alignment that carries a small deletion straddling
// an exon boundary lifts to an N whose position is off by the deletion length: the boundary sits on a
// non-canonical motif a few bp from the true junction (e.g. AGRN 15M140N where 18M140N is the GT-AG
// placement). Re-partitioning the flanking M|N|M moves matched read bases across the intron without
// changing read or reference span, so the fix is a pure CIGAR edit - validated against the read so a
// slide is only taken when every moved base still matches the reference (no new mismatch) and the
// motif strictly improves. The alignment start never moves (only interior boundaries shift).
//
// Runs after tail extension as the last lifted-cigar pass, so it operates on the final boundaries.
// Gated on a reference being available; junctions already on a canonical motif are skipped cheaply.
public class JunctionCanonicalizer
{
    // max bp an intron is slid to reach a canonical motif. Bounds the search to the size of a tx-contig
    // boundary deletion (SPLICE_FLANKING_DELETION_MAX_BP), the only place the lifted junction drifts.
    public static final int DEFAULT_MAX_SHIFT = 5;

    private final RefSequenceSource mRefSource;
    private final int mMaxShift;
    private final AtomicLong mJunctionsShifted = new AtomicLong();

    public JunctionCanonicalizer(final RefSequenceSource refSource, final int maxShift)
    {
        mRefSource = refSource;
        mMaxShift = maxShift;
    }

    public long junctionsShifted() { return mJunctionsShifted.get(); }

    public JunctionCanonicalizationResult tryCanonicalize(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
            return JunctionCanonicalizationResult.unchanged();

        final List<CigarShape.Element> elements = CigarShape.parse(cigar);
        if(elements.size() < 3)
            return JunctionCanonicalizationResult.unchanged();

        int shifted = 0;

        // walk left-to-right; for each interior N with matched flanks, attempt a canonicalizing slide.
        // Read/genome anchors are recomputed from the live list each iteration so an applied shift (which
        // re-sizes the flanking M's) is reflected when the next junction is evaluated.
        for(int j = 1; j < elements.size() - 1; ++j)
        {
            if(elements.get(j).Op != OP_SKIPPED)
                continue;
            if(!isMatchedOp(elements.get(j - 1).Op) || !isMatchedOp(elements.get(j + 1).Op))
                continue;

            if(trySlideJunction(chromosome, alignmentStart, elements, j, readBases))
                ++shifted;
        }

        if(shifted == 0)
            return JunctionCanonicalizationResult.unchanged();

        mJunctionsShifted.addAndGet(shifted);
        return new JunctionCanonicalizationResult(true, CigarShape.format(elements), shifted);
    }

    // Evaluates intron at element index j and applies the best canonicalizing shift in place. Returns
    // true if a shift was applied.
    private boolean trySlideJunction(
            final String chromosome, final int alignmentStart, final List<CigarShape.Element> elements,
            final int j, final byte[] readBases)
    {
        // genome position of the first base after the left flank = intron start; read index of the
        // first base of the right flank, computed by summing live element lengths up to j.
        int genomePos = alignmentStart;
        int readIndex = 0;
        for(int i = 0; i < j; ++i)
        {
            final CigarShape.Element e = elements.get(i);
            if(consumesReference(e.Op))
                genomePos += e.Length;
            if(consumesRead(e.Op))
                readIndex += e.Length;
        }

        final int intronLength = elements.get(j).Length;
        final int intronStart = genomePos;                 // first intronic base (1-based)
        final int intronEnd = genomePos + intronLength - 1; // last intronic base
        final int rightReadStart = readIndex;              // read index of right flank's first base

        final int leftLen = elements.get(j - 1).Length;
        final int rightLen = elements.get(j + 1).Length;
        final int maxShift = Math.min(mMaxShift, Math.min(leftLen - 1, rightLen - 1));
        if(maxShift < 1)
            return false;

        // ref windows wide enough to read the donor/acceptor motif and the moved bases at any shift.
        final int donorBlockStart = intronStart - mMaxShift;
        final byte[] donorBlock = mRefSource.getBases(chromosome, donorBlockStart, intronStart + mMaxShift + 1);
        final int acceptorBlockStart = intronEnd - mMaxShift - 1;
        final byte[] acceptorBlock = mRefSource.getBases(chromosome, acceptorBlockStart, intronEnd + mMaxShift);
        if(donorBlock == null || acceptorBlock == null)
            return false;

        final int currentTier = motifTier(donorBlock, donorBlockStart, intronStart, acceptorBlock, acceptorBlockStart, intronEnd, 0);
        if(currentTier >= SpliceMotif.TIER_CANONICAL)
            return false;

        int bestShift = 0;
        int bestTier = currentTier;
        // smallest |shift| wins ties; try in increasing magnitude, positive before negative.
        for(int mag = 1; mag <= maxShift; ++mag)
        {
            for(final int d : new int[] {mag, -mag})
            {
                final int tier = motifTier(
                        donorBlock, donorBlockStart, intronStart, acceptorBlock, acceptorBlockStart, intronEnd, d);
                if(tier <= bestTier)
                    continue;
                if(!movedBasesMatch(d, intronStart, intronEnd, rightReadStart, readBases, donorBlock, donorBlockStart,
                        acceptorBlock, acceptorBlockStart))
                    continue;
                bestTier = tier;
                bestShift = d;
            }
            if(bestShift != 0)
                break; // smallest magnitude that improves the motif is preferred
        }

        if(bestShift == 0)
            return false;

        // apply: +d grows the left flank and shrinks the right; -d the reverse. Intron length unchanged.
        elements.set(j - 1, new CigarShape.Element(leftLen + bestShift, elements.get(j - 1).Op));
        elements.set(j + 1, new CigarShape.Element(rightLen - bestShift, elements.get(j + 1).Op));
        return true;
    }

    // splice-motif tier of the intron when shifted by d: donor at [intronStart+d, +1], acceptor at
    // [intronEnd+d-1, intronEnd+d].
    private static int motifTier(
            final byte[] donorBlock, final int donorBlockStart, final int intronStart,
            final byte[] acceptorBlock, final int acceptorBlockStart, final int intronEnd, final int d)
    {
        final byte[] donor = twoBases(donorBlock, intronStart + d - donorBlockStart);
        final byte[] acceptor = twoBases(acceptorBlock, intronEnd + d - 1 - acceptorBlockStart);
        return SpliceMotif.classify(donor, acceptor);
    }

    // verifies the bases that cross the intron under shift d still match the reference at their new
    // positions. +d moves the right flank's first d bases onto the donor side (genome [intronStart,
    // intronStart+d-1]); -d moves the left flank's last d bases onto the acceptor side (genome
    // [intronEnd+d+1, intronEnd]). Requires an exact match for every moved base.
    private static boolean movedBasesMatch(
            final int d, final int intronStart, final int intronEnd, final int rightReadStart, final byte[] readBases,
            final byte[] donorBlock, final int donorBlockStart, final byte[] acceptorBlock, final int acceptorBlockStart)
    {
        if(d > 0)
        {
            for(int k = 0; k < d; ++k)
            {
                final int readIdx = rightReadStart + k;
                if(readIdx >= readBases.length)
                    return false;
                final int refIdx = intronStart + k - donorBlockStart;
                if(refIdx < 0 || refIdx >= donorBlock.length)
                    return false;
                if(!basesEqualIgnoreCase(readBases[readIdx], donorBlock[refIdx]))
                    return false;
            }
            return true;
        }

        final int shift = -d;
        for(int k = 0; k < shift; ++k)
        {
            final int readIdx = rightReadStart - shift + k;
            if(readIdx < 0)
                return false;
            final int genome = intronEnd + d + 1 + k; // = intronEnd - shift + 1 + k
            final int refIdx = genome - acceptorBlockStart;
            if(refIdx < 0 || refIdx >= acceptorBlock.length)
                return false;
            if(!basesEqualIgnoreCase(readBases[readIdx], acceptorBlock[refIdx]))
                return false;
        }
        return true;
    }

    private static byte[] twoBases(final byte[] block, final int offset)
    {
        if(offset < 0 || offset + 1 >= block.length)
            return null;
        return new byte[] {block[offset], block[offset + 1]};
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
}
