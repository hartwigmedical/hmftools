package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class RegionData extends ChrBaseRegion
{
    public final String GeneName;
    public final RegionType Type;
    public final int ExonRank;

    private String mExtraInfo;

    public RegionData(final String geneName, final ChrBaseRegion region, final int exonRank, final RegionType type)
    {
        super(region.Chromosome, region.start(), region.end());
        GeneName = geneName;
        Type = type;
        ExonRank = exonRank;
        mExtraInfo = "";
    }

    public void setExtraInfo(final String extraInfo) { mExtraInfo = extraInfo; }
    public String getExtraInfo() { return mExtraInfo; }

    public String name()
    {
        if(Type == RegionType.CODING)
            return String.format("%s_%d_%s", GeneName, ExonRank, Type);

        if(Type == RegionType.INTRONIC)
            return String.format("%s_%s_%s", GeneName, Type, mExtraInfo);

        int posMidpoint = (start() + end()) / 2;
        return String.format("%s_%s_%d", GeneName, Chromosome, posMidpoint);
    }

    public String idName() { return name(); }

    public String toString()
    {
        return String.format("%s: %s:%d_%d %s", GeneName, Chromosome, start(), end(), Type);
    }

    public static RegionData fromSpecificRegionCsv(final String data)
    {
        final String[] items = data.split(",", -1);

        // Chromosome,PosStart,PosEnd,GeneName,Type,Info
        ChrBaseRegion region = new ChrBaseRegion(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]));
        RegionData regionData = new RegionData(items[3], region, 0, RegionType.valueOf(items[4]));
        regionData.setExtraInfo(items[5]);
        return regionData;
    }

    public static void mergeRegion(final List<RegionData> regions, final RegionData newRegion)
    {
        // insert regions in ascending order by position
        // merge any overlapping regions
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(newRegion.start() > region.end())
            {
                ++index;
                continue;
            }

            if(region.start() > newRegion.end())
                break;

            if(newRegion.matches(region))
                return;

            // handle merges
            int startPosition = min(region.start(), newRegion.start());
            region.setStart(startPosition);

            int endPosition = max(region.end(), newRegion.end());
            region.setEnd(endPosition);

            ++index;

            while(index < regions.size())
            {
                RegionData nextRegion = regions.get(index);

                if(nextRegion.start() > region.end())
                    break;

                endPosition = max(region.end(), nextRegion.end());
                region.setEnd(endPosition);
                regions.remove(index);
            }

            return;
        }

        regions.add(index, newRegion);
    }

    public static void integrateRegion(final List<RegionData> regions, final RegionData newRegion)
    {
        // split these new regions if they overlap an existing coding region
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(newRegion.start() > region.end())
            {
                ++index;
                continue;
            }

            if(region.start() > newRegion.end())
                break;

            if(newRegion.matches(region))
                return;

            if(newRegion.start() < region.start())
            {
                RegionData preRegion = new RegionData(
                        newRegion.GeneName,
                        new ChrBaseRegion(newRegion.Chromosome, newRegion.start(), region.start() - 1),
                        newRegion.ExonRank, newRegion.Type);

                preRegion.setExtraInfo(newRegion.getExtraInfo());

                regions.add(index, preRegion);
                ++index; // for the additional insert

                // adjust for remaining segment
                if(newRegion.end() <= region.end())
                    return;

                newRegion.setStart(region.end() + 1);
            }
            else
            {
                newRegion.setStart(region.end() + 1);
            }

            ++index;
        }

        regions.add(index, newRegion);
    }

    public static boolean validate(final List<RegionData> regions)
    {
        for(int i = 0; i < regions.size() - 1; ++i)
        {
            RegionData region = regions.get(i);
            RegionData nextRegion = regions.get(i + 1);

            if(nextRegion.start() <= region.end())
                return false;
        }

        return true;
    }
}
