package com.hartwig.hmftools.esvee.caller.annotation;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion.SPARE_FIELD;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PonSglRegion implements Comparable<PonSglRegion>
{
    public final ChrBaseRegion Region;
    public final Orientation Orient;
    public final int PonCount;

    public PonSglRegion(final ChrBaseRegion region, final Orientation orient, final int ponCount)
    {
        Region = region;
        Orient = orient;
        PonCount = ponCount;
    }

    public boolean overlaps(final BaseRegion svStart)
    {
        return positionsOverlap(Region.start(), Region.end(), svStart.start(), svStart.end());
    }

    public boolean matches(final BaseRegion svRegion, Orientation orientation)
    {
        return overlaps(svRegion) && Orient == orientation;
    }

    public String toString()
    {
        return String.format("region(%s) orient(%d) pon(%d)", Region, Orient.asByte(), PonCount);
    }

    @Override
    public int compareTo(final PonSglRegion other)
    {
        if(Region.start() == other.Region.start())
        {
            if(Region.end() == other.Region.end())
                return 0;

            return Region.end() < other.Region.end() ? -1 : 1;
        }

        return Region.start() < other.Region.start() ? -1 : 1;
    }

    public static PonSglRegion fromBedRecord(final String data)
    {
        String[] items = data.split(TSV_DELIM, -1);

        ChrBaseRegion region = new ChrBaseRegion(items[0], Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));

        Orientation orient = Orientation.fromChar(items[5].charAt(0));
        int ponCount = Integer.parseInt(items[4]);

        return new PonSglRegion(region, orient, ponCount);
    }

    public String toBedRecord()
    {
        // re-applying the 1 start offset

        // fields: Chr,PosBegin,PosEnd,Unknown,PonCount,Orientation
        return String.format("%s\t%d\t%d\t%s\t%d\t%s",
                Region.Chromosome, Region.start() - 1, Region.end(), SPARE_FIELD, PonCount, Orient.asChar());

    }
}
