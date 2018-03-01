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

public class LineElementAnnotator {

    private List<GenomeRegion> mLineElements;
    private static final int PERMITTED_DISTANCE = 2000;

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public LineElementAnnotator()
    {
        mLineElements = Lists.newArrayList();
    }

    public void loadLineElementsFile(final String filename)
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

                mLineElements.add(genomeRegion);

                LOGGER.debug("loaded line element: chr({}) pos({}-{})",
                        genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
            }

            LOGGER.debug("loaded {} line element records", mLineElements.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read line element CSV file({})", mLineElements);
        }
    }

    public boolean isLineElement(final SvClusterData svData, final boolean useStart)
    {
        if(mLineElements.isEmpty())
            return false;

        for(final GenomeRegion genomeRegion : mLineElements)
        {
            if(!genomeRegion.chromosome().equals(svData.chromosome(useStart)))
                continue;

            // test if the SV falls within the LE +/- a buffer
            if(svData.position(useStart) >= genomeRegion.start() - PERMITTED_DISTANCE
            && svData.position(useStart) <= genomeRegion.end() + PERMITTED_DISTANCE)
            {
                LOGGER.debug("var({}) found in line element({}.>{})",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return true;
            }
        }

        return false;
    }

}
