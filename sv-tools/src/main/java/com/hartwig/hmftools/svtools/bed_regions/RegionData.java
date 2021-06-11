package com.hartwig.hmftools.svtools.bed_regions;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.List;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class RegionData
{
    public final String GeneName;
    public final BaseRegion Region;
    public final RegionType Type;
    private String mExtraInfo;

    public RegionData(final String geneName, final BaseRegion region, final RegionType type)
    {
        GeneName = geneName;
        Region = region;
        Type = type;
        mExtraInfo = "";
    }

    public void setExtraInfo(final String extraInfo) { mExtraInfo = extraInfo; }
    public String getExtraInfo() { return mExtraInfo; }

    public String name()
    {
        if(Type == RegionType.CODING)
            return String.format("%s_%s", GeneName, Type);

        if(Type == RegionType.INTRONIC)
            return String.format("%s_%s_%s", GeneName, Type, mExtraInfo);

        int posMidpoint = (Region.start() + Region.end()) / 2;
        return String.format("%s_%s_%d", GeneName, Region.Chromosome, posMidpoint);
    }

    public String toString()
    {
        return String.format("%s: %s %s", GeneName, Region, Type);
    }

    public static RegionData fromSpecificRegionCsv(final String data)
    {
        final String[] items = data.split(",", -1   );

        // Chromosome,PosStart,PosEnd,GeneName,Type,Info
        BaseRegion region = new BaseRegion(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]));
        RegionData regionData = new RegionData(items[3], region, RegionType.valueOf(items[4]));
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

            if(newRegion.Region.start() > region.Region.end())
            {
                ++index;
                continue;
            }

            if(region.Region.start() > newRegion.Region.end())
                break;

            if(newRegion.Region.matches(region.Region))
                return;

            // handle merges
            int startPosition = min(region.Region.start(), newRegion.Region.start());
            region.Region.setStart(startPosition);

            int endPosition = max(region.Region.end(), newRegion.Region.end());
            region.Region.setEnd(endPosition);

            ++index;

            while(index < regions.size())
            {
                RegionData nextRegion = regions.get(index);

                if(nextRegion.Region.start() > region.Region.end())
                    break;

                endPosition = max(region.Region.end(), nextRegion.Region.end());
                region.Region.setEnd(endPosition);
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

            if(newRegion.Region.start() > region.Region.end())
            {
                ++index;
                continue;
            }

            if(region.Region.start() > newRegion.Region.end())
                break;

            if(newRegion.Region.matches(region.Region))
                return;

            if(newRegion.Region.start() < region.Region.start())
            {
                RegionData preRegion = new RegionData(
                        newRegion.GeneName,
                        new BaseRegion(newRegion.Region.Chromosome, newRegion.Region.start(), region.Region.start() - 1),
                        newRegion.Type);

                preRegion.setExtraInfo(newRegion.getExtraInfo());

                regions.add(index, preRegion);
                ++index; // for the additional insert

                // adjust for remaining segment
                if(newRegion.Region.end() <= region.Region.end())
                    return;

                newRegion.Region.setStart(region.Region.end() + 1);
            }
            else
            {
                newRegion.Region.setStart(region.Region.end() + 1);
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

            if(nextRegion.Region.start() <= region.Region.end())
                return false;
        }

        return true;
    }
}
