package com.hartwig.hmftools.tars.liftback.features;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

import htsjdk.samtools.SAMRecord;

// Recomputes bwa-style alignment scores after candidates are lifted to genome space. Cigars are left unchanged: the only
// soft-clip walk in TARS is the weak-overhang collapse in OverhangGate.
public class GenomicAlignmentScorer
{
    private final RefGenomeInterface mRefGenome;

    public GenomicAlignmentScorer(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public boolean enabled()
    {
        return mRefGenome != null;
    }

    public void scoreCandidates(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(!enabled() || alignments.size() < 2)
        {
            return;
        }

        byte[] forwardBases = record.getReadBases();
        if(forwardBases == null || forwardBases.length == 0)
        {
            return;
        }

        boolean recordForward = !record.getReadNegativeStrandFlag();
        byte[] reverseBases = hasOppositeStrandCandidate(alignments, recordForward) ? reverseComplement(forwardBases) : null;

        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.Dropped || alignment.LiftedCigar == null)
            {
                continue;
            }

            score(alignment, basesFor(alignment, recordForward, forwardBases, reverseBases));
        }
    }

    public void score(final LiftedAlignment alignment, final byte[] readBases)
    {
        if(alignment == null || alignment.Dropped || alignment.LiftedCigar == null)
        {
            return;
        }

        alignment.GenomicScore = BwaScoring.genomicScore(
                mRefGenome, alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, readBases);
    }

    public int scoreRecord(final LiftedAlignment alignment, final SAMRecord record)
    {
        if(!enabled() || alignment == null || alignment.Dropped || alignment.LiftedCigar == null)
        {
            return Integer.MIN_VALUE;
        }

        byte[] bases = record.getReadBases();
        if(bases == null || bases.length == 0)
        {
            return Integer.MIN_VALUE;
        }
        if(alignment.ForwardStrand == record.getReadNegativeStrandFlag())
        {
            bases = reverseComplement(bases);
        }
        return BwaScoring.genomicScore(
                mRefGenome, alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, bases);
    }

    private static boolean hasOppositeStrandCandidate(final List<LiftedAlignment> alignments, final boolean recordForward)
    {
        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.ForwardStrand != recordForward)
            {
                return true;
            }
        }
        return false;
    }

    private static byte[] basesFor(
            final LiftedAlignment alignment, final boolean recordForward, final byte[] forwardBases,
            final byte[] reverseBases)
    {
        return alignment.ForwardStrand == recordForward ? forwardBases : reverseBases;
    }

    private static byte[] reverseComplement(final byte[] bases)
    {
        byte[] reversed = Arrays.copyOf(bases, bases.length);
        Nucleotides.reverseComplementBasesInPlace(reversed, 0, reversed.length);
        return reversed;
    }
}
