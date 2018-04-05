package com.hartwig.hmftools.svanalysis.annotators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FragileSiteAnnotator {

    private List<GenomeRegion> mFragileSites;
    private List<GenomeRegion> mIdentifiedFragileSites;

    public static final String KNOWN_FS = "true";
    public static final String IDENTIFIED_FS = "ident";
    public static final String NO_FS = "false";

    private static final String CSV_FS_TYPE_IDENTIFIED = "Identified";

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public FragileSiteAnnotator()
    {
        mFragileSites = Lists.newArrayList();
        mIdentifiedFragileSites = Lists.newArrayList();
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

                if(items.length < 7)
                    continue;

                final GenomeRegion genomeRegion = GenomeRegionFactory.create(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]));

                if(items[6].equals(CSV_FS_TYPE_IDENTIFIED))
                    mIdentifiedFragileSites.add(genomeRegion);
                else
                    mFragileSites.add(genomeRegion);

//                LOGGER.debug("loaded fragile site: chr({}) pos({}-{})",
//                        genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
            }

            LOGGER.debug("loaded {} known fragile site records, {} identifed", mFragileSites.size(), mIdentifiedFragileSites.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read fragile site CSV file({})", filename);
        }
    }

    public String isFragileSite(final SvClusterData svData, final boolean useStart)
    {
        if(mFragileSites.isEmpty())
            return NO_FS;

        for(final GenomeRegion genomeRegion : mFragileSites)
        {
            if(genomeRegion.chromosome().equals(svData.chromosome(useStart))
            && genomeRegion.start() <= svData.position(useStart) && genomeRegion.end() >= svData.position(useStart))
            {
                LOGGER.debug("var({}) found in known fragile site",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return KNOWN_FS;
            }
        }

        for(final GenomeRegion genomeRegion : mIdentifiedFragileSites)
        {
            if(genomeRegion.chromosome().equals(svData.chromosome(useStart))
            && genomeRegion.start() <= svData.position(useStart) && genomeRegion.end() >= svData.position(useStart))
            {
                LOGGER.debug("var({}) found in identified fragile site",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return IDENTIFIED_FS;
            }
        }

        return NO_FS;
    }

}
