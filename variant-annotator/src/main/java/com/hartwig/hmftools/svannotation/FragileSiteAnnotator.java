package com.hartwig.hmftools.svannotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.svannotation.analysis.SvClusterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FragileSiteAnnotator {

    private List<GenomeRegion> mFragileSites;

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public FragileSiteAnnotator()
    {
        mFragileSites = Lists.newArrayList();
    }

    public void loadFragileSitesFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try {

            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line;
            while ((line = fileReader.readLine()) != null) {

                if(line.contains("Chromosome"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < 3)
                    continue;

                final GenomeRegion genomeRegion = GenomeRegionFactory.create(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]));

                mFragileSites.add(genomeRegion);

                LOGGER.debug("loaded fragile site: chr({}) pos({}-{})",
                        genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
            }

            LOGGER.debug("loaded {} fragile site records", mFragileSites.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read fragile site CSV file({})", mFragileSites);
        }
    }

    public boolean isFragileSite(final SvClusterData svData, final boolean useStart)
    {
        if(mFragileSites.isEmpty())
            return false;

        for(final GenomeRegion genomeRegion : mFragileSites)
        {
            if(genomeRegion.chromosome().equals(svData.chromosome(useStart))
            && genomeRegion.start() <= svData.position(useStart) && genomeRegion.end() >= svData.position(useStart))
            {
                LOGGER.debug("var({}) found in fragile site",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return true;
            }
        }

        return false;
    }

}
