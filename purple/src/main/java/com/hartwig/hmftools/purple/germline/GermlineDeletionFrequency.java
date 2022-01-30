package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.purple.drivers.DeletionRegionFrequency;

public class GermlineDeletionFrequency
{
    private Map<String,List<DeletionRegionFrequency>> mChrRegionMap;

    public static final String COHORT_DEL_FREQ_FILE = "germline_del_freq_file";

    public GermlineDeletionFrequency(final String filename)
    {
        mChrRegionMap = Maps.newHashMap();
        loadCohortFrequencies(filename);
    }

    public int getRegionFrequency(final String chromsome, int regionStart, int regionEnd, int buffer)
    {
        List<DeletionRegionFrequency> regions = mChrRegionMap.get(chromsome);

        if(regions == null)
            return 0;

        int totalFrequency = 0;

        for(DeletionRegionFrequency regionFreq : regions)
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
            int frequencyIndex = fieldsIndexMap.get("Frequency");

            String currentChromosome = "";
            List<DeletionRegionFrequency> regions = null;

            for(String line : lines)
            {
                final String[] values = line.split(",", -1);

                String chromosome = values[chromosomeIndex];
                int regionStart = Integer.parseInt(values[regionStartIndex]);
                int regionEnd = Integer.parseInt(values[regionEndIndex]);
                int frequency = Integer.parseInt(values[frequencyIndex]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    regions = Lists.newArrayList();
                    mChrRegionMap.put(chromosome, regions);
                }

                regions.add(new DeletionRegionFrequency(new BaseRegion(regionStart, regionEnd), frequency));
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
