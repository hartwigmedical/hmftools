package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.FLD_FREQUENCY;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.FLD_TYPE;

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
    private Map<String, List<AmpDelRegionFrequency>> mChrRegionMap;

    public static final String COHORT_AMP_DEL_FREQ_FILE = "`germline_amp_del_freq_file`";

    public GermlineAmpDelFrequencyCache(final String filename)
    {
        mChrRegionMap = Maps.newHashMap();
        loadCohortFrequencies(filename);
    }

    public int getRegionFrequency(final String chromsome, int regionStart, int regionEnd, int buffer,
            AmpDelRegionFrequency.EventType eventType)
    {
        if(!mChrRegionMap.containsKey(chromsome))
        {
            return 0;
        }
        List<AmpDelRegionFrequency> regions = mChrRegionMap.get(chromsome).stream().filter(r -> r.Type == eventType).toList();

        int totalFrequency = 0;

        for(AmpDelRegionFrequency regionFreq : regions)
        {
            if(regionStart - buffer > regionFreq.Region.start())
            {
                continue;
            }

            if(regionFreq.Region.start() - buffer > regionStart)
            {
                break;
            }

            if(abs(regionFreq.Region.start() - regionStart) <= buffer && abs(regionFreq.Region.end() - regionEnd) <= buffer)
            {
                totalFrequency += regionFreq.Frequency;
            }
        }

        return totalFrequency;
    }

    private void loadCohortFrequencies(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);

            int chromosomeIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int regionStartIndex = fieldsIndexMap.get(FLD_REGION_START);
            int regionEndIndex = fieldsIndexMap.get(FLD_REGION_END);
            int typeIndex = fieldsIndexMap.getOrDefault(FLD_TYPE, -1);
            int frequencyIndex = fieldsIndexMap.get(FLD_FREQUENCY);

            String currentChromosome = "";
            List<AmpDelRegionFrequency> regions = null;

            for(String line : lines)
            {
                final String[] values = line.split(CSV_DELIM, -1);

                String chromosome = values[chromosomeIndex];
                int regionStart = Integer.parseInt(values[regionStartIndex]);
                int regionEnd = Integer.parseInt(values[regionEndIndex]);

                AmpDelRegionFrequency.EventType type = typeIndex >= 0
                        ? AmpDelRegionFrequency.EventType.valueOf(values[typeIndex]) : AmpDelRegionFrequency.EventType.DEL;

                int frequency = Integer.parseInt(values[frequencyIndex]);

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;

                    regions = mChrRegionMap.get(chromosome);

                    if(regions == null)
                    {
                        regions = Lists.newArrayList();
                        mChrRegionMap.put(chromosome, regions);
                    }
                }

                regions.add(new AmpDelRegionFrequency(new BaseRegion(regionStart, regionEnd), type, frequency));
            }

            PPL_LOGGER.info("loaded {} germline amp-del frequencies from file: {}", lines.size(), filename);
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to read cohort germline amp-del file: {}", e.toString());
        }
    }
}
