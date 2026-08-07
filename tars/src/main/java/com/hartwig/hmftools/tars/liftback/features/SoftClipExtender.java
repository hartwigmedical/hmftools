package com.hartwig.hmftools.tars.liftback.features;

import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.tars.common.BwaScoring.maxScoringPrefix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.common.BwaScoring;
import com.hartwig.hmftools.tars.common.TarsCigarUtils;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SoftClipExtender
{
    private final RefGenomeInterface mRefGenome;

    public SoftClipExtender(final RefGenomeInterface refGenome)
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

    public LiftedAlignment extend(final LiftedAlignment alignment, final byte[] readBases)
    {
        if(alignment == null || alignment.Dropped || alignment.LiftedCigar == null)
        {
            return alignment;
        }

        Result extended = extend(alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar, readBases);
        if(extended.pos() == alignment.LiftedPos && extended.cigar().equals(alignment.LiftedCigar))
        {
            return alignment;
        }
        return alignment.withLiftedCigar(extended.pos(), extended.cigar());
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

    private Result extend(final String chromosome, final int alignmentStart, final String cigar, final byte[] readBases)
    {
        if(!enabled() || chromosome == null || cigar == null || readBases == null || readBases.length == 0)
        {
            return new Result(alignmentStart, cigar);
        }

        List<CigarElement> elements = CigarUtils.cigarElementsFromStr(cigar);
        if(elements.isEmpty() || CigarUtils.hasHardClip(elements))
        {
            return new Result(alignmentStart, cigar);
        }

        List<CigarElement> working = new ArrayList<>(elements);
        int alignmentEnd = alignmentStart + CigarUtils.cigarAlignedLength(working) - 1;

        int trailing = extendSide(false, chromosome, alignmentStart, alignmentEnd, readBases, working);
        int leading = extendSide(true, chromosome, alignmentStart, alignmentEnd, readBases, working);

        if(trailing == 0 && leading == 0)
        {
            return new Result(alignmentStart, cigar);
        }

        return new Result(alignmentStart - leading, CigarUtils.cigarElementsToStr(working));
    }

    private int extendSide(
            final boolean leading, final String chromosome, final int alignmentStart, final int alignmentEnd,
            final byte[] readBases, final List<CigarElement> working)
    {
        if(working.size() < 2)
        {
            return 0;
        }

        int softclipIndex = leading ? 0 : working.size() - 1;
        int matchedIndex = leading ? 1 : working.size() - 2;

        CigarElement softclip = working.get(softclipIndex);
        if(softclip.getOperator() != CigarOperator.S || !TarsCigarUtils.isMatchedOp(working.get(matchedIndex).getOperator()))
        {
            return 0;
        }

        int budget = leading ? Math.min(softclip.getLength(), alignmentStart - 1) : softclip.getLength();
        if(budget <= 0)
        {
            return 0;
        }

        byte[] refBases = leading
                ? refBases(chromosome, alignmentStart - budget, alignmentStart - 1)
                : refBases(chromosome, alignmentEnd + 1, alignmentEnd + budget);
        if(refBases == null || refBases.length == 0)
        {
            return 0;
        }

        int walkLength = Math.min(refBases.length, budget);
        int softclipLength = softclip.getLength();
        int extended = leading
                ? maxScoringPrefix(
                        reverseArray(Arrays.copyOfRange(readBases, softclipLength - walkLength, softclipLength)),
                        reverseArray(Arrays.copyOfRange(refBases, refBases.length - walkLength, refBases.length)))
                : maxScoringPrefix(
                        Arrays.copyOfRange(readBases, readBases.length - softclipLength, readBases.length - softclipLength + walkLength),
                        Arrays.copyOfRange(refBases, 0, walkLength));

        if(extended == 0)
        {
            return 0;
        }

        TarsCigarUtils.extendTerminalSoftClipIntoMatch(working, leading, extended);
        return extended;
    }

    private byte[] refBases(final String chromosome, final int posStart, final int posEnd)
    {
        return BwaScoring.refWindow(mRefGenome, chromosome, posStart, posEnd);
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

    private record Result(int pos, String cigar)
    {
    }
}
