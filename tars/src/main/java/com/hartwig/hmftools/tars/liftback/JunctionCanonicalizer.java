package com.hartwig.hmftools.tars.liftback;

import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.SpliceMotif;
import com.hartwig.hmftools.tars.liftback.rescue.Tier;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

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
    private final RefSequenceSource mRefSource;
    private final int mMaxShift;
    private long mJunctionsShifted;

    public JunctionCanonicalizer(final RefSequenceSource refSource, final int maxShift)
    {
        mRefSource = refSource;
        mMaxShift = maxShift;
    }

    public long junctionsShifted() { return mJunctionsShifted; }

    // Restore the counter to roll back a discarded provisional mate decision without double-counting
    // (see LiftBackGroupProcessor.processNameGroup).
    public void restoreJunctionsShifted(final long junctionsShifted) { mJunctionsShifted = junctionsShifted; }

    public JunctionCanonicalizationResult tryCanonicalize(
            final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(mRefSource == null || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return JunctionCanonicalizationResult.unchanged();
        }

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.size() < 3)
        {
            return JunctionCanonicalizationResult.unchanged();
        }

        int shifted = 0;

        // walk left-to-right; for each interior N with matched flanks, attempt a canonicalizing slide.
        // Read/genome anchors are recomputed from the live list each iteration so an applied shift (which
        // re-sizes the flanking M's) is reflected when the next junction is evaluated.
        for(int j = 1; j < elements.size() - 1; ++j)
        {
            if(elements.get(j).getOperator() != CigarOperator.N)
                continue;
            if(!isMatchedOp(elements.get(j - 1).getOperator()) || !isMatchedOp(elements.get(j + 1).getOperator()))
                continue;

            if(trySlideJunction(chromosome, alignmentStart, elements, j, readBases))
            {
                ++shifted;
            }
        }

        if(shifted == 0)
        {
            return JunctionCanonicalizationResult.unchanged();
        }

        mJunctionsShifted += shifted;
        return new JunctionCanonicalizationResult(true, CigarUtils.cigarElementsToStr(elements));
    }

    // Evaluates intron at element index j and applies the best canonicalizing shift in place. Returns
    // true if a shift was applied.
    private boolean trySlideJunction(
            final String chromosome, final int alignmentStart, final List<CigarElement> elements,
            final int j, final byte[] readBases)
    {
        // genome position of the first base after the left flank = intron start; read index of the
        // first base of the right flank, computed by summing live element lengths up to j.
        int genomePos = alignmentStart;
        int readIndex = 0;
        for(int i = 0; i < j; ++i)
        {
            CigarElement e = elements.get(i);
            if(e.getOperator().consumesReferenceBases())
            {
                genomePos += e.getLength();
            }
            if(e.getOperator().consumesReadBases())
            {
                readIndex += e.getLength();
            }
        }

        int intronLength = elements.get(j).getLength();
        int intronStart = genomePos;                 // first intronic base (1-based)
        int intronEnd = genomePos + intronLength - 1; // last intronic base
        int rightReadStart = readIndex;              // read index of right flank's first base

        int leftLen = elements.get(j - 1).getLength();
        int rightLen = elements.get(j + 1).getLength();
        int maxShift = Math.min(mMaxShift, Math.min(leftLen - 1, rightLen - 1));
        if(maxShift < 1)
        {
            return false;
        }

        // ref windows wide enough to read the donor/acceptor motif and the moved bases at any shift.
        int donorBlockStart = intronStart - mMaxShift;
        byte[] donorBlock = mRefSource.getBases(chromosome, donorBlockStart, intronStart + mMaxShift + 1);
        int acceptorBlockStart = intronEnd - mMaxShift - 1;
        byte[] acceptorBlock = mRefSource.getBases(chromosome, acceptorBlockStart, intronEnd + mMaxShift);
        if(donorBlock == null || acceptorBlock == null)
        {
            return false;
        }

        Tier currentTier = motifTier(donorBlock, donorBlockStart, intronStart, acceptorBlock, acceptorBlockStart, intronEnd, 0);
        if(currentTier.compareTo(Tier.CANONICAL) >= 0)
        {
            return false;
        }

        int bestShift = 0;
        Tier bestTier = currentTier;
        // smallest |shift| wins ties; try in increasing magnitude, positive before negative.
        for(int mag = 1; mag <= maxShift; ++mag)
        {
            for(final int d : new int[] { mag, -mag })
            {
                Tier tier = motifTier(
                        donorBlock, donorBlockStart, intronStart, acceptorBlock, acceptorBlockStart, intronEnd, d);
                if(tier.compareTo(bestTier) <= 0)
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
        {
            return false;
        }

        // apply: +d grows the left flank and shrinks the right; -d the reverse. Intron length unchanged.
        elements.set(j - 1, new CigarElement(leftLen + bestShift, elements.get(j - 1).getOperator()));
        elements.set(j + 1, new CigarElement(rightLen - bestShift, elements.get(j + 1).getOperator()));
        return true;
    }

    // splice-motif tier of the intron when shifted by d: donor at [intronStart+d, +1], acceptor at
    // [intronEnd+d-1, intronEnd+d].
    private static Tier motifTier(
            final byte[] donorBlock, final int donorBlockStart, final int intronStart,
            final byte[] acceptorBlock, final int acceptorBlockStart, final int intronEnd, final int d)
    {
        byte[] donor = twoBases(donorBlock, intronStart + d - donorBlockStart);
        byte[] acceptor = twoBases(acceptorBlock, intronEnd + d - 1 - acceptorBlockStart);
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
                int readIdx = rightReadStart + k;
                if(readIdx >= readBases.length)
                {
                    return false;
                }
                int refIdx = intronStart + k - donorBlockStart;
                if(refIdx < 0 || refIdx >= donorBlock.length)
                {
                    return false;
                }
                if(!basesEqualIgnoreCase(readBases[readIdx], donorBlock[refIdx]))
                {
                    return false;
                }
            }
            return true;
        }

        int shift = -d;
        for(int k = 0; k < shift; ++k)
        {
            int readIdx = rightReadStart - shift + k;
            if(readIdx < 0)
            {
                return false;
            }
            int genome = intronEnd + d + 1 + k; // = intronEnd - shift + 1 + k
            int refIdx = genome - acceptorBlockStart;
            if(refIdx < 0 || refIdx >= acceptorBlock.length)
            {
                return false;
            }
            if(!basesEqualIgnoreCase(readBases[readIdx], acceptorBlock[refIdx]))
            {
                return false;
            }
        }
        return true;
    }

    private static byte[] twoBases(final byte[] block, final int offset)
    {
        if(offset < 0 || offset + 1 >= block.length)
        {
            return null;
        }
        return new byte[] { block[offset], block[offset + 1] };
    }

    private static boolean isMatchedOp(final CigarOperator op)
    {
        return op == CigarOperator.M || op == CigarOperator.EQ || op == CigarOperator.X;
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        if(a == b)
        {
            return true;
        }
        return (a & ~0x20) == (b & ~0x20);
    }
}
