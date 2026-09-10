package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.tars.common.TarsCigarUtils.mergeAdjacentSameOp;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

// Cigars are left unchanged: the only soft-clip walk in TARS is the weak-overhang collapse in OverhangGate.
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

    public void scorePlacements(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(!enabled() || alignments.size() < 2)
        {
            return;
        }

        ReadBases readBases = ReadBases.of(record);
        if(readBases == null)
        {
            return;
        }

        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.Dropped || alignment.LiftedCigar == null)
            {
                continue;
            }

            score(alignment, readBases.forAlignment(alignment));
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

        ReadBases readBases = ReadBases.of(record);
        if(readBases == null)
        {
            return Integer.MIN_VALUE;
        }

        return BwaScoring.genomicScore(
                mRefGenome, alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar,
                readBases.forAlignment(alignment));
    }

    // An intron-only lift preserves BWA's scoring operations: the transcript contig contains the same exonic bases and
    // N has no alignment-score cost. Any other CIGAR change makes the input AS stale and requires a genomic rescore.
    public static boolean requiresRescore(final SAMRecord record, final LiftedAlignment alignment)
    {
        List<CigarElement> liftedOps = new ArrayList<>();
        for(CigarElement element : TextCigarCodec.decode(alignment.LiftedCigar).getCigarElements())
        {
            if(element.getOperator() != CigarOperator.N)
            {
                liftedOps.add(element);
            }
        }

        List<CigarElement> inputOps = new ArrayList<>(record.getCigar().getCigarElements());
        if(record.getReadNegativeStrandFlag() == alignment.ForwardStrand)
        {
            Collections.reverse(inputOps);
        }

        return !mergeAdjacentSameOp(liftedOps).equals(mergeAdjacentSameOp(inputOps));
    }
}
