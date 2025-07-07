package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.Double.NaN;
import static java.lang.String.format;

import java.util.Optional;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PanelRegion extends ChrBaseRegion
{
    public final RegionType Type;
    public final String SourceInfo;

    public final String Sequence;
    public final double GcContent;
    public final Optional<Double> QualityScore;

    public final boolean IsProbe;

    public PanelRegion(
            final ChrBaseRegion region, final RegionType type, final String sourceInfo,
            final String sequence, final double gcContent, final double qualityScore)
    {
        super(region.Chromosome, region.start(), region.end());
        Type = type;
        SourceInfo = sourceInfo;

        Sequence = sequence;
        GcContent = gcContent;
        QualityScore = Optional.of(qualityScore);
        IsProbe = true;
    }

    public PanelRegion(final ChrBaseRegion region, final RegionType type, final String sourceInfo)
    {
        super(region.Chromosome, region.start(), region.end());
        Type = type;
        SourceInfo = sourceInfo;

        Sequence = "";
        GcContent = -1;
        QualityScore = Optional.empty();
        IsProbe = false;
    }

    public String extendedSourceInfo()
    {
        if(!IsProbe)
        {
            return SourceInfo;
        }

        return format("%s;gc=%.3f,qualityScore=%.0f", SourceInfo, GcContent, QualityScore.orElse(NaN));
    }

    public String toString()
    {
        if(IsProbe)
        {
            return format("region(%s) type(%s) info(%s)", super.toString(), Type, SourceInfo);
        }
        else
        {
            return format("region(%s) type(%s) info(%s) qualityScore(%.3f)", super.toString(), Type, SourceInfo, QualityScore.orElse(NaN));
        }
    }
}
