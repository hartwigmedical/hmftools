package com.hartwig.hmftools.esvee.caller.annotation;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PonSvRegion implements Comparable<PonSvRegion>
{
    public final ChrBaseRegion RegionStart;
    public final Orientation OrientStart;
    public final ChrBaseRegion RegionEnd;
    public final Orientation OrientEnd;
    public final int PonCount;

    protected static final String SPARE_FIELD = ".";

    public PonSvRegion(
            final ChrBaseRegion regionStart, final Orientation orientStart,
            final ChrBaseRegion regionEnd, final Orientation orientEnd, final int ponCount)
    {
        RegionStart = regionStart;
        OrientStart = orientStart;
        RegionEnd = regionEnd;
        OrientEnd = orientEnd;
        PonCount = ponCount;
    }

    public boolean overlapsStart(final BaseRegion svStart)
    {
        return positionsOverlap(RegionStart.start(), RegionStart.end(), svStart.start(), svStart.end());
    }

    public boolean matches(final BaseRegion svStart, final ChrBaseRegion svEnd, Orientation orientStart, Orientation orientEnd)
    {
        return overlapsStart(svStart) && RegionEnd.overlaps(svEnd) && OrientStart == orientStart && OrientEnd == orientEnd;
    }

    @Override
    public int compareTo(final PonSvRegion other)
    {
        if(RegionStart.start() == other.RegionStart.start())
        {
            if(RegionStart.end() == other.RegionStart.end())
                return 0;

            return RegionStart.end() < other.RegionStart.end() ? -1 : 1;
        }

        return RegionStart.start() < other.RegionStart.start() ? -1 : 1;
    }

    public static PonSvRegion fromBedRecord(final String data)
    {
        String[] items = data.split(TSV_DELIM, -1);

        ChrBaseRegion regionStart = new ChrBaseRegion(items[0], Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));
        ChrBaseRegion regionEnd = new ChrBaseRegion(items[3], Integer.parseInt(items[4]) + 1, Integer.parseInt(items[5]));

        Orientation orientStart = Orientation.fromChar(items[8].charAt(0));
        Orientation orientEnd = Orientation.fromChar(items[9].charAt(0));
        int ponCount = Integer.parseInt(items[7]);

        return new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, ponCount);
    }

    public String toBedRecord()
    {
        // re-applying the 1 start offset
        // fields: ChrStart,PosStartBegin,PosStartEnd,ChrEnd,PosEndBegin,PosEndEnd,Unknown,PonCount,OrientStart,OrientEnd
        return String.format("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s",
                RegionStart.Chromosome, RegionStart.start() - 1, RegionStart.end(),
                RegionEnd.Chromosome, RegionEnd.start() - 1, RegionEnd.end(), SPARE_FIELD,
                PonCount, OrientStart.asChar(), OrientEnd.asChar());
    }

    public String toString()
    {
        return String.format("start(%s:%d) end(%s:%d) pon(%d)",
                RegionStart, OrientStart.asByte(), RegionEnd, OrientEnd.asByte(), PonCount);
    }
}
