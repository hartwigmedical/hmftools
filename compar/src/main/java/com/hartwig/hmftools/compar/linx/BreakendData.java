package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class BreakendData
{
    public final LinxBreakend Breakend;

    // SV info
    public final String SvId;
    public final StructuralVariantType SvType;
    public final String Chromosome;
    public final int Position;
    public final int[] HomologyOffset;
    public final byte Orientation;
    public final String ComparisonChromosome;
    public final int ComparisonPosition;

    public BreakendData(
            final LinxBreakend breakend, final String svId, final StructuralVariantType svType, final String chromosome,
            final int position, final byte orientation, final int[] homologyOffset, final String comparisonChromosome,
            final int comparisonPosition)
    {
        Breakend = breakend;
        SvId = svId;
        SvType = svType;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        HomologyOffset = homologyOffset;
        ComparisonChromosome = comparisonChromosome;
        ComparisonPosition = comparisonPosition;
    }

    public boolean matches(final BreakendData other)
    {
        // check SV info
        if(SvType != other.SvType || !ComparisonChromosome.equals(other.Chromosome) || Orientation != other.Orientation)
            return false;

        // check breakend info
        if(!positionsOverlap(
                ComparisonPosition + HomologyOffset[0], ComparisonPosition + HomologyOffset[1],
                other.Position + other.HomologyOffset[0], other.Position + other.HomologyOffset[1]))
        {
            return false;
        }

        boolean checkTranscriptIds = !Breakend.canonical() || !other.Breakend.canonical();
        if(checkTranscriptIds && !other.Breakend.transcriptId().equals(Breakend.transcriptId()))
        {
            return false;
        }

        return true;
    }

    public String svInfoStr()
    {
        if(!ComparisonChromosome.equals(Chromosome) || ComparisonPosition != Position)
        {
            return format("%s:%s %s:%d:%d (%s:%d)", SvId, SvType, Chromosome, Position, Orientation, ComparisonChromosome, ComparisonPosition);
        }
        else
        {
            return format("%s:%s %s:%d:%d", SvId, SvType, Chromosome, Position, Orientation);
        }
    }

    public String transcriptStr()
    {
        return format("%s:%s:%s:%d", Breakend.transcriptId(), Breakend.codingType(), Breakend.regionType(), Breakend.exonUp());
    }

    public String fullStr()
    {
        return format("%s reported(%s) transcript(%s)", svInfoStr(), Breakend.reportedStatus(), transcriptStr());
    }
}
