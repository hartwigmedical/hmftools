package com.hartwig.hmftools.isofox.unmapped;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class UnmappedRead
{
    public final String ReadId;
    public final ChrBaseRegion ReadRegion;
    public final byte Orientation;
    public final int ScLength;
    public final int ScSide;
    public final double AvgBaseQual;
    public final String GeneName;
    public final String TransName;
    public final int ExonRank;
    public final int ExonBoundary;
    public final int ExonDistance;
    public final String SpliceType;
    public final String ScBases;

    public UnmappedRead(
            final String readId, final ChrBaseRegion readRegion, final byte orientation, final int scLength, final int scSide,
            final double avgBaseQual, final String geneName, final String transName, final int exonRank, final int exonBoundary,
            final int exonDistance, final String spliceType, final String scBases)
    {
        ReadId = readId;
        ReadRegion = readRegion;
        Orientation = orientation;
        ScLength = scLength;
        ScSide = scSide;
        AvgBaseQual = avgBaseQual;
        GeneName = geneName;
        TransName = transName;
        ExonRank = exonRank;
        ExonBoundary = exonBoundary;
        ExonDistance = exonDistance;
        SpliceType = spliceType;
        ScBases = scBases;
    }

    public String formKey()
    {
        return String.format("%s_%d_%d_%d", ReadRegion.Chromosome, ExonBoundary, ScSide, Orientation);
    }

    public boolean matches(final UnmappedRead other)
    {
        return ReadRegion.Chromosome.equals(other.ReadRegion.Chromosome) && ExonBoundary == other.ExonBoundary
                && ScSide == other.ScSide && Orientation == other.Orientation;
    }
}
