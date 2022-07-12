package com.hartwig.hmftools.svprep;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

public class BlacklistLocations
{
    private final Map<String,List<BaseRegion>> mChrLocationsMap; // keyed by chromosome start

    public BlacklistLocations(final String filename)
    {
        mChrLocationsMap = Maps.newHashMap();
        loadFile(filename);
    }

    public List<BaseRegion> getRegions(final String chromosome) { return mChrLocationsMap.get(chromosome); }

    public boolean inBlacklistLocation(final String chromosome, final int position)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(position);
        return regions != null ? regions.stream().anyMatch(x -> x.containsPosition(position)) : false;
    }

    public void addRegion(final String chromosome, final BaseRegion region)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrLocationsMap.put(chromosome, regions);
        }

        regions.add(region);
    }

    private void sortAndMerge()
    {
        for(List<BaseRegion> regions : mChrLocationsMap.values())
        {
            Collections.sort(regions);

            // merge any adjacent regions
            int index = 0;
            while(index < regions.size() - 1)
            {
                BaseRegion region = regions.get(index);
                BaseRegion nextRegion = regions.get(index + 1);

                if(region.end() >= nextRegion.start() - 2)
                {
                    region.setEnd(nextRegion.end());
                    regions.remove(index + 1);
                }
                else
                {
                    ++index;
                }
            }
        }
    }

    private void loadFile(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            int itemCount = 0;
            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                if(line.contains("chromosome") && line.contains("regionType"))
                    continue;

                final String[] values = line.split("\t", -1);

                if(values.length < 3)
                {
                    SV_LOGGER.error("invalid blacklist entry: {}", line);
                    return;
                }

                String chromosome = values[0];
                int posStart = Integer.parseInt(values[1]) + 1;
                int posEnd = Integer.parseInt(values[2]);

                addRegion(chromosome, new BaseRegion(posStart, posEnd));
                ++itemCount;
            }

            SV_LOGGER.info("loaded {} blacklist locations from file", itemCount, filename);

            sortAndMerge();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to load hotspot data file({}): {}", filename, e.toString());
            return;
        }
    }
}
