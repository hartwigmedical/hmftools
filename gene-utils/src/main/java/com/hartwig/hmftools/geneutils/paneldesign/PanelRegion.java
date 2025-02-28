package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PanelRegion extends ChrBaseRegion
{
    public final RegionType Type;
    public final String SourceInfo;

    public final String Sequence;
    public final double GcContent;
    public final double BlastnScore;

    public final boolean IsProbe;

    public PanelRegion(
            final ChrBaseRegion region, final RegionType type, final String sourceInfo,
            final String sequence, final double gcContent, final double blastnScore)
    {
        super(region.Chromosome, region.start(), region.end());
        Type = type;
        SourceInfo = sourceInfo;

        Sequence = sequence;
        GcContent = gcContent;
        BlastnScore = blastnScore;
        IsProbe = true;
    }

    public PanelRegion(final ChrBaseRegion region, final RegionType type, final String sourceInfo)
    {
        super(region.Chromosome, region.start(), region.end());
        Type = type;
        SourceInfo = sourceInfo;

        Sequence = "";
        GcContent = -1;
        BlastnScore = -1;
        IsProbe = false;
    }

    public String extendedSourceInfo()
    {
        if(!IsProbe)
            return SourceInfo;

        return format("%s;gc=%.3f,blastScore=%.0f", SourceInfo, GcContent, BlastnScore);
    }

    public String toString()
    {
        if(IsProbe)
            return format("region(%s) type(%s) info(%s)", super.toString(), Type, SourceInfo);
        else
            return format("region(%s) type(%s) info(%s) blastScore(%.3f)", super.toString(), Type, SourceInfo, BlastnScore);
    }
}
