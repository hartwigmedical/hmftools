package com.hartwig.hmftools.geneutils.fusion;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.lowerChromosome;

public class FusionBedData implements Comparable<FusionBedData>
{
    public final String Name;
    public final String ChrUp;
    public final String ChrDown;
    public final byte StrandUp;
    public final byte StrandDown;
    public final int PositionUpStart;
    public final int PositionUpEnd;
    public final int PositionDownStart;
    public final int PositionDownEnd;

    private final boolean mUpIsStart;

    public FusionBedData(
            final String name, final String chrUp, final String chrDown, final byte strandUp, final byte strandDown,
            final int positionUpStart, final int positionUpEnd, final int positionDownStart, final int positionDownEnd)
    {
        Name = name;
        ChrUp = chrUp;
        ChrDown = chrDown;
        StrandUp = strandUp;
        StrandDown = strandDown;
        PositionUpStart = positionUpStart;
        PositionUpEnd = positionUpEnd;
        PositionDownStart = positionDownStart;
        PositionDownEnd = positionDownEnd;

        if(lowerChromosome(ChrUp, ChrDown))
        {
            mUpIsStart = true;
        }
        else if(ChrUp.equals(ChrDown))
        {
            mUpIsStart = adjustedUpStart() < adjustedDownStart();
        }
        else
        {
            mUpIsStart = false;
        }
    }

    // methods for sorting
    public int adjustedUpStart()
    {
        return StrandUp == POS_STRAND ? PositionUpStart - GenerateFusionFiles.PRE_GENE_BUFFER - 1 : PositionUpStart - 1;
    }

    public int adjustedUpEnd()
    {
        return StrandUp == POS_STRAND ? PositionUpEnd : PositionUpEnd + GenerateFusionFiles.PRE_GENE_BUFFER;
    }

    public int adjustedDownStart()
    {
        return StrandDown == POS_STRAND ? PositionDownStart - GenerateFusionFiles.PRE_GENE_BUFFER - 1 : PositionDownStart - 1;
    }

    public int adjustedDownEnd()
    {
        return StrandDown == POS_STRAND ? PositionDownEnd : PositionDownEnd + GenerateFusionFiles.PRE_GENE_BUFFER;
    }

    public char strandUpChar()
    {
        return StrandUp == POS_STRAND ? '+' : '-';
    }

    // REVERSE STRAND2 since for the downstream genes the orientation is opposite to upstream (ie +ve strand = -ve orientation and vice versa)
    public char strandDownChar()
    {
        return StrandDown == POS_STRAND ? '-' : '+';
    }

    // ordered for the BED
    public String chrStart()
    {
        return mUpIsStart ? ChrUp : ChrDown;
    }

    public String chrEnd()
    {
        return !mUpIsStart ? ChrUp : ChrDown;
    }

    public int posStartStart()
    {
        return mUpIsStart ? adjustedUpStart() : adjustedDownStart();
    }

    public int posStartEnd()
    {
        return mUpIsStart ? adjustedUpEnd() : adjustedDownEnd();
    }

    public int posEndStart()
    {
        return !mUpIsStart ? adjustedUpStart() : adjustedDownStart();
    }

    public int posEndEnd()
    {
        return !mUpIsStart ? adjustedUpEnd() : adjustedDownEnd();
    }

    public char strandStart()
    {
        return mUpIsStart ? strandUpChar() : strandDownChar();
    }

    public char strandEnd()
    {
        return !mUpIsStart ? strandUpChar() : strandDownChar();
    }

    public String toString()
    {
        return format("%s start(%s:%d-%d) end(%s:%d-%d) upIsStart(%s)",
                Name, chrStart(), posStartStart(), posStartEnd(), chrEnd(), posEndStart(), posEndEnd(), mUpIsStart);
    }

    @Override
    public int compareTo(final FusionBedData other)
    {
        if(lowerChromosome(chrStart(), other.chrStart()))
        {
            return -1;
        }
        else if(chrStart().equals(other.chrStart()))
        {
            if(posStartStart() == other.posStartStart())
            {
                return Name.compareTo(other.Name);
            }

            return posStartStart() < other.posStartStart() ? -1 : 1;
        }
        else
        {
            return 1;
        }
    }
}
