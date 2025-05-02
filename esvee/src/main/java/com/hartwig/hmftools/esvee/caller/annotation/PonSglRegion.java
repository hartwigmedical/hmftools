package com.hartwig.hmftools.esvee.caller.annotation;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.esvee.caller.annotation.PonCache.FLD_PON_COUNT;
import static com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion.SPARE_FIELD;

import java.util.StringJoiner;

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
        if(Region.start() != other.Region.start())
            return Region.start() < other.Region.start() ? -1 : 1;

        if(Region.end() != other.Region.end())
            return Region.end() < other.Region.end() ? -1 : 1;

        if(Orient != other.Orient)
            return Orient.isForward() ? -1 : 1;

        return Integer.compare(PonCount, other.PonCount);
    }

    public static PonSglRegion fromTsv(final String data)
    {
        String[] items = data.split(TSV_DELIM, -1);

        int index = 0;
        ChrBaseRegion region = new ChrBaseRegion(items[index++], Integer.parseInt(items[index++]), Integer.parseInt(items[index++]));

        Orientation orientation = Orientation.fromByteStr(items[index++]);
        int ponCount = Integer.parseInt(items[index]);

        return new PonSglRegion(region, orientation, ponCount);
    }

    public static PonSglRegion fromBedRecord(final String data)
    {
        String[] items = data.split(TSV_DELIM, -1);

        ChrBaseRegion region = new ChrBaseRegion(items[0], Integer.parseInt(items[1]) + 1, Integer.parseInt(items[2]));

        Orientation orient = Orientation.fromChar(items[5].charAt(0));
        int ponCount = Integer.parseInt(items[4]);

        return new PonSglRegion(region, orient, ponCount);
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END).add(FLD_ORIENTATION).add(FLD_PON_COUNT);
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(Region.Chromosome).add(String.valueOf(Region.start())).add(String.valueOf(Region.end()));
        sj.add(String.valueOf(Orient.asByte())).add(String.valueOf(PonCount));
        return sj.toString();
    }

    public String toBedRecord()
    {
        // re-applying the 1 start offset

        // fields: Chr,PosBegin,PosEnd,Unknown,PonCount,Orientation
        return String.format("%s\t%d\t%d\t%s\t%d\t%s",
                Region.Chromosome, Region.start() - 1, Region.end(), SPARE_FIELD, PonCount, Orient.asChar());

    }
}
