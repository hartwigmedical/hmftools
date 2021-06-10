package com.hartwig.hmftools.svtools.bed_regions;

import static java.lang.Math.max;

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

    public void setExtraInfo(final String extraInfo)
    {
        mExtraInfo = extraInfo;
    }

    public String name()
    {
        if(Type == RegionType.CODING)
            return String.format("%s_%s", GeneName, Type);
        else
            return String.format("%s_%s_%s", GeneName, Type, mExtraInfo);
    }

    public String toString()
    {
        return String.format("%s: %s %s", GeneName, Region, Type);
    }

    public static RegionData fromSpecificRegionCsv(final String data)
    {
        final String[] items = data.split(",");

        // Chromosome,PosStart,PosEnd,GeneName,Type,Info
        BaseRegion region = new BaseRegion(items[0], Integer.parseInt(items[1]), Integer.parseInt(items[2]));
        RegionData regionData = new RegionData(items[3], region, RegionType.valueOf(items[4]));
        regionData.setExtraInfo(items[5]);
        return regionData;
    }

    public static void addRegion(final List<RegionData> regions, final RegionData newRegionData)
    {
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(newRegionData.Region.start() > region.Region.end())
            {
                ++index;
                continue;
            }

            if(region.Region.start() > newRegionData.Region.end())
                break;

            // handle merges
            int startPosition = max(region.Region.start(), newRegionData.Region.start());
            region.Region.setStart(startPosition);

            int endPosition = max(region.Region.end(), newRegionData.Region.end());
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

        regions.add(index, newRegionData);
    }

    public static boolean validate(final List<RegionData> regions)
    {
        for(int i = 0; i < regions.size(); ++i)
        {
            RegionData region1 = regions.get(i);

            for(int j = i + 1; j < regions.size(); ++j)
            {
                RegionData region2 = regions.get(j);

                if(region1.Region.overlaps(region2.Region))
                    return false;
            }
        }

        return true;
    }
}
