package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.Orientation;

public class FragmentData
{
    public final PrepRead Read;
    public final int MatePosition;
    public final boolean IsPrimary;
    public final boolean PrimaryIsLower;

    public FragmentData(final PrepRead read)
    {
        Read = read;

        IsPrimary = !read.isSupplementaryAlignment();

        SupplementaryReadData suppData = read.supplementaryAlignment();

        if(read.Chromosome.equals(suppData.Chromosome))
        {
            PrimaryIsLower = read.start() <= suppData.Position;
        }
        else
        {
            PrimaryIsLower = HumanChromosome.lowerChromosome(read.Chromosome, suppData.Chromosome);
        }

        if(read.record().hasAttribute(MATE_CIGAR_ATTRIBUTE))
        {
            String mateCigar = read.record().getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            MatePosition = getFivePrimeUnclippedPosition(read.record().getMateAlignmentStart(), mateCigar, Read.mateOrientation().isForward());
        }
        else
        {
            // fall-back, imprecise duplicate evaluation
            MatePosition = read.MatePosStart;
        }
    }

    public String mateChromosome() { return Read.MateChromosome; }
    public Orientation mateOrientation() { return Read.mateOrientation(); }

    public String toString()
    {
        return format("id(%s) mate(%s:%d:%d) %s %s",
                Read.id(), mateChromosome(), MatePosition, mateOrientation().asByte(),
                IsPrimary ? "primary" : "supp", PrimaryIsLower ? "primary-lower" : "primary-upper");
    }
}
