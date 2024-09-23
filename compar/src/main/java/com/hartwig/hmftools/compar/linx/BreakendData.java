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

    public BreakendData(
            final LinxBreakend breakend, final String svId, final StructuralVariantType svType, final String chromosome,
            final int position, final byte orientation, final int[] homologyOffset)
    {
        Breakend = breakend;
        SvId = svId;
        SvType = svType;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        HomologyOffset = homologyOffset;
    }

    public boolean matches(final BreakendData other)
    {
        // check SV info
        if(SvType != other.SvType || !Chromosome.equals(other.Chromosome) || Orientation != other.Orientation)
            return false;

        // check breakend info
        if(!positionsOverlap(
                Position + HomologyOffset[0], Position + HomologyOffset[1],
                other.Position + other.HomologyOffset[0], other.Position + other.HomologyOffset[1]))
        {
            return false;
        }

        return Breakend.transcriptId().equals(other.Breakend.transcriptId());
    }

    public String svInfoStr()
    {
        return format("%s:%s %s:%d:%d", SvId, SvType, Chromosome, Position, Orientation);
    }

    public String transcriptStr()
    {
        return format("%s:%s:%d", Breakend.codingType(), Breakend.regionType(), Breakend.exonUp());
    }

    public String fullStr()
    {
        return format("%s reported(%d) transcript(%s)", svInfoStr(), Breakend.reportedDisruption(), transcriptStr());
    }
}
