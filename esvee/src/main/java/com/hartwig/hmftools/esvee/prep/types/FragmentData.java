package com.hartwig.hmftools.esvee.prep.types;

import static java.lang.Math.abs;
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
    public final boolean IsPrimary;
    public final boolean ReadIsLowerVsSuppData;

    public final int OtherPosition; // mate unclipped 5' position, or 3' end for unpaired reads

    public FragmentData(final PrepRead read)
    {
        Read = read;

        IsPrimary = !read.isSupplementaryAlignment();

        SupplementaryReadData suppData = read.supplementaryAlignment();

        if(read.Chromosome.equals(suppData.Chromosome))
        {
            ReadIsLowerVsSuppData = read.AlignmentStart <= suppData.Position;
        }
        else
        {
            ReadIsLowerVsSuppData = HumanChromosome.lowerChromosome(read.Chromosome, suppData.Chromosome);
        }

        if(read.isPaired())
        {
            if(read.record().hasAttribute(MATE_CIGAR_ATTRIBUTE))
            {
                String mateCigar = read.record().getStringAttribute(MATE_CIGAR_ATTRIBUTE);
                OtherPosition = getFivePrimeUnclippedPosition(
                        read.record().getMateAlignmentStart(), mateCigar, Read.mateOrientation().isForward());
            }
            else
            {
                // fall-back, imprecise duplicate evaluation
                OtherPosition = read.MatePosStart;
            }
        }
        else
        {
            OtherPosition = read.orientation().isForward() ? read.UnclippedEnd : read.UnclippedStart;
        }
    }

    public String otherChromosome() { return Read.isPaired() ? Read.MateChromosome : Read.Chromosome; }
    public Orientation otherOrientation() { return Read.isPaired() ? Read.mateOrientation() : Read.orientation(); }

    public boolean isConsensus() { return Read.record().hasAttribute(CONSENSUS_READ_ATTRIBUTE); }

    public boolean matches(final FragmentData other)
    {
        return otherChromosome().equals(other.otherChromosome())
                && otherOrientation() == other.otherOrientation()
                && OtherPosition == other.OtherPosition;
    }

    public boolean withinPositionRange(final FragmentData other, final int permittedPositionDiff)
    {
        if(Read.orientation() != other.Read.orientation())
            return false;

        int startDiff = abs(Read.UnclippedStart - other.Read.UnclippedStart);
        int endDiff = abs(Read.UnclippedEnd - other.Read.UnclippedEnd);
        return startDiff + endDiff <= permittedPositionDiff;
    }

    public static int unclippedPosition(final PrepRead read)
    {
        return read.orientation().isForward() ? read.UnclippedStart : read.UnclippedEnd;
    }

    public String toString()
    {
        return format("id(%s) mate(%s:%d:%d) %s %s",
                Read.id(), otherChromosome(), OtherPosition, otherOrientation().asByte(),
                IsPrimary ? "primary" : "supp", ReadIsLowerVsSuppData ? "lower" : "upper");
    }
}
