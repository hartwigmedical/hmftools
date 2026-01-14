package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;

public class GermlineAmpDelFrequencyCache
{
    private Map<String,List<AmpDelRegionFrequency>> mChrRegionMap;

    public static final String COHORT_AMP_DEL_FREQ_FILE = "germline_amp_del_freq_file";

    public GermlineAmpDelFrequencyCache(final String filename)
    {
        mChrRegionMap = Maps.newHashMap();
        loadCohortFrequencies(filename);
    }

    public int getRegionFrequency(final String chromsome, int regionStart, int regionEnd, int buffer, AmpDelRegionFrequency.EventType eventType)
    {
        if (!mChrRegionMap.containsKey(chromsome))
        {
            return 0;
        }
        List<AmpDelRegionFrequency> regions = mChrRegionMap.get(chromsome).stream().filter(r -> r.Type == eventType).toList();

        int totalFrequency = 0;

        for(AmpDelRegionFrequency regionFreq : regions)
        {
            if(regionStart - buffer > regionFreq.Region.start())
                continue;

            if(regionFreq.Region.start() - buffer > regionStart)
                break;

            if(abs(regionFreq.Region.start() - regionStart) <= buffer && abs(regionFreq.Region.end() - regionEnd) <= buffer)
                totalFrequency += regionFreq.Frequency;
        }

        return totalFrequency;
    }

    private boolean loadCohortFrequencies(final String filename)
    {
        if(filename == null)
            return false;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");

            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int regionStartIndex = fieldsIndexMap.get("RegionStart");
            int regionEndIndex = fieldsIndexMap.get("RegionEnd");
            int typeIndex = fieldsIndexMap.get("Type");
            int frequencyIndex = fieldsIndexMap.get("Frequency");

            String currentChromosome = "";
            List<AmpDelRegionFrequency> regions = null;

            for(String line : lines)
            {
                final String[] values = line.split(",", -1);

                String chromosome = values[chromosomeIndex];
                int regionStart = Integer.parseInt(values[regionStartIndex]);
                int regionEnd = Integer.parseInt(values[regionEndIndex]);
                AmpDelRegionFrequency.EventType type = AmpDelRegionFrequency.EventType.valueOf(values[typeIndex]);
                int frequency = Integer.parseInt(values[frequencyIndex]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    regions = Lists.newArrayList();
                    mChrRegionMap.put(chromosome, regions);
                }

                regions.add(new AmpDelRegionFrequency(new BaseRegion(regionStart, regionEnd), type, frequency));
            }

            PPL_LOGGER.info("loaded {} germline deletions frequencies from file: {}", lines.size(), filename);
            return true;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to read cohort germline deletions file: {}", e.toString());
            return false;
        }
    }
}
