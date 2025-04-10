package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNC_POS_DIFF;

import com.hartwig.hmftools.esvee.assembly.types.Junction;

public final class ReadFilters
{
    public static boolean recordSoftClipsAndCrossesJunction(final Read read, final Junction junction)
    {
        int junctionOverlap = 0;

        if(junction.isForward())
        {
            // first check indel-inferred soft-clips
            if(read.hasIndelImpliedUnclippedEnd())
                return read.maxUnclippedEnd() > junction.Position;

            // soft-clip must be close enough to the junction
            if(read.isRightClipped())
                return abs(read.alignmentEnd() - junction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF && read.unclippedEnd() > junction.Position;

            // consider alignments with small junction overlap, which ought to have mismatched bases
            junctionOverlap = read.alignmentEnd() - junction.Position;
        }
        else
        {
            if(read.hasIndelImpliedUnclippedStart())
                return read.minUnclippedStart() < junction.Position;

            if(read.isLeftClipped())
                return abs(read.alignmentStart() - junction.Position) <= ASSEMBLY_MAX_JUNC_POS_DIFF && read.unclippedStart() < junction.Position;

            junctionOverlap = junction.Position - read.alignmentStart();
        }

        return read.snvCount() > 0 && junctionOverlap > 0 && junctionOverlap <= ASSEMBLY_MAX_JUNC_POS_DIFF;
    }

    public static boolean recordSoftClipsAtJunction(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() && read.alignmentEnd() == junction.Position;
        }
        else
        {
            return read.isLeftClipped() && read.alignmentStart() == junction.Position;
        }
    }

    public static int readJunctionExtensionLength(final Read read, final Junction junction)
    {
        if(junction.isForward())
        {
            return read.isRightClipped() ? max(read.maxUnclippedEnd() - junction.Position, 0) : 0;
        }
        else
        {
            return read.isLeftClipped() ? max(junction.Position - read.minUnclippedStart(), 0) : 0;
        }
    }
}
