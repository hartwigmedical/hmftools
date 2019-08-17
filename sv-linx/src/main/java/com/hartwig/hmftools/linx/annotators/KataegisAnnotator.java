package com.hartwig.hmftools.linx.annotators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class KataegisAnnotator
{
    private Map<String,Map<String,List<KataegisData>>> mSampleChrData;

    private static int CSV_REQUIRED_FIELDS = 6;

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public KataegisAnnotator()
    {
        mSampleChrData = Maps.newHashMap();
    }

    public void loadKataegisData(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int recordCount = 0;

            // SampleId,KataegisId,SnvCount,Chromosome,MinPosition,MaxPosition

            String line;
            String currentSample = "";
            Map<String,List<KataegisData>> chrData = null;

            while ((line = fileReader.readLine()) != null)
            {
                if(line.contains("SampleId"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < CSV_REQUIRED_FIELDS)
                    continue;

                final String sampleId = items[0];
                final String chromosome = items[3];

                KataegisData data = new KataegisData(
                        chromosome,
                        Long.parseLong(items[4]),
                        Long.parseLong(items[5]),
                        items[1],
                        Integer.parseInt(items[2]));

                if(!currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    chrData = Maps.newHashMap();
                    mSampleChrData.put(sampleId, chrData);
                }

                List<KataegisData> dataList = chrData.get(chromosome);

                if(dataList == null)
                {
                    dataList = Lists.newArrayList();
                    chrData.put(chromosome, dataList);
                }

                dataList.add(data);
                ++recordCount;
            }

            LOGGER.debug("loaded {} kataegis data records, samples({})", recordCount, mSampleChrData.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read kataegis CSV file({})", filename);
        }
    }

    /*
    public boolean isFragileSite(final SvVarData svData, final boolean useStart)
    {
        if(mFragileSites.isEmpty())
            return false;

        for(final GenomeRegion genomeRegion : mFragileSites)
        {
            if(genomeRegion.chromosome().equals(svData.chromosome(useStart))
                    && genomeRegion.start() <= svData.position(useStart) && genomeRegion.end() >= svData.position(useStart))
            {
                LOGGER.debug("var({}) found in known fragile site",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return true;
            }
        }

        return false;
    }
    */

}
