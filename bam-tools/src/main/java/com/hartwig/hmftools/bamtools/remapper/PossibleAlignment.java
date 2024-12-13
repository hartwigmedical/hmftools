package com.hartwig.hmftools.bamtools.remapper;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.esvee.assembly.alignment.AlternativeAlignment;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;

public class PossibleAlignment implements Comparable<PossibleAlignment>
{
    private final BwaMemAlignment mateAlignment;
    private final AlternativeAlignment alternativeAlignment;

    public PossibleAlignment(final BwaMemAlignment mateAlignment, final AlternativeAlignment alternativeAlignment)
    {
        this.mateAlignment = mateAlignment;
        this.alternativeAlignment = alternativeAlignment;
        if(mateAlignment.getRefId() != 5)
        {
            throw new IllegalArgumentException("possible alignments only considered for chr6");
        }
        if(HumanChromosome.fromString(alternativeAlignment.Chromosome) != HumanChromosome._6)
        {
            throw new IllegalArgumentException("possible alignments only considered for chr6");
        }
        if(!areDirectionallyCompatible(mateAlignment, alternativeAlignment))
        {
            throw new IllegalArgumentException(
                    "Not directionally compatible: " + mateAlignment.getSamFlag() + " and " + alternativeAlignment);
        }
    }

    static boolean areDirectionallyCompatible(BwaMemAlignment alignment, AlternativeAlignment alternativeAlignment)
    {
        if(SAMFlag.getFlags(alignment.getSamFlag()).contains(SAMFlag.READ_REVERSE_STRAND))
        {
            return alternativeAlignment.Orient.isForward();
        }
        return alternativeAlignment.Orient.isReverse();
    }

    public int distanceFromMate()
    {
        return Math.abs(Math.abs(mateAlignment.getRefStart()) - Math.abs(alternativeAlignment.Position));
    }

    public AlternativeAlignment alignment()
    {
        return alternativeAlignment;
    }

    @Override
    public int compareTo(@NotNull final PossibleAlignment o)
    {
        return distanceFromMate() - o.distanceFromMate();
    }
}
