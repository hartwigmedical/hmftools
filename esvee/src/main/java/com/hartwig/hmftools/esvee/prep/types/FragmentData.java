package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
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
    public final boolean ReadIsLowerVsSuppData;

    public FragmentData(final PrepRead read)
    {
        Read = read;

        IsPrimary = !read.isSupplementaryAlignment();

        SupplementaryReadData suppData = read.supplementaryAlignment();

        if(read.Chromosome.equals(suppData.Chromosome))
        {
            ReadIsLowerVsSuppData = read.start() <= suppData.Position;
        }
        else
        {
            ReadIsLowerVsSuppData = HumanChromosome.lowerChromosome(read.Chromosome, suppData.Chromosome);
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

    public boolean isConsensus() { return Read.record().hasAttribute(CONSENSUS_READ_ATTRIBUTE); }

    public boolean matches(final FragmentData other)
    {
        return mateChromosome().equals(other.mateChromosome())
                && mateOrientation() == other.mateOrientation()
                && MatePosition == other.MatePosition;
    }

    public static int unclippedPosition(final PrepRead read)
    {
        return read.orientation().isForward() ? read.unclippedStart() : read.unclippedEnd();
    }

    public String toString()
    {
        return format("id(%s) mate(%s:%d:%d) %s %s",
                Read.id(), mateChromosome(), MatePosition, mateOrientation().asByte(),
                IsPrimary ? "primary" : "supp", ReadIsLowerVsSuppData ? "lower" : "upper");
    }
}
