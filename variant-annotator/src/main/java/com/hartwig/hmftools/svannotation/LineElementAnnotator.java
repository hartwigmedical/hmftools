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

    public static final String KNOWN_LINE_ELEMENT = "true";
    public static final String IDENTIFIED_LINE_ELEMENT = "ident";
    public static final String NO_LINE_ELEMENT = "false";

    private static final String CSV_LE_TYPE_IDENTIFIED = "Identified";

    private List<GenomeRegion> mKnownLineElements;
    private List<GenomeRegion> mIdentifiedLineElements;
    private static final int PERMITTED_DISTANCE = 5000;

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public LineElementAnnotator()
    {
        mKnownLineElements = Lists.newArrayList();
        mIdentifiedLineElements = Lists.newArrayList();
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

                if(items.length < 4)
                    continue;

                final GenomeRegion genomeRegion = GenomeRegionFactory.create(items[0], Long.parseLong(items[1]), Long.parseLong(items[2]));

                if(items[3].equals(CSV_LE_TYPE_IDENTIFIED))
                    mIdentifiedLineElements.add(genomeRegion);
                else
                    mKnownLineElements.add(genomeRegion);

//                LOGGER.debug("loaded line element: chr({}) pos({}-{})",
//                        genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
            }

            LOGGER.debug("loaded {} known line elements, {} identified", mKnownLineElements.size(), mIdentifiedLineElements.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read line element CSV file({})", filename);
        }
    }

    public String isLineElement(final SvClusterData svData, final boolean useStart)
    {
        if(mKnownLineElements.isEmpty() && mIdentifiedLineElements.isEmpty())
            return NO_LINE_ELEMENT;

        for(final GenomeRegion genomeRegion : mKnownLineElements)
        {
            if(!genomeRegion.chromosome().equals(svData.chromosome(useStart)))
                continue;

            // test if the SV falls within the LE +/- a buffer
            if(svData.position(useStart) >= genomeRegion.start() - PERMITTED_DISTANCE
            && svData.position(useStart) <= genomeRegion.end() + PERMITTED_DISTANCE)
            {
                LOGGER.debug("var({}) found in known line element({}.>{})",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return KNOWN_LINE_ELEMENT;
            }
        }

        for(final GenomeRegion genomeRegion : mIdentifiedLineElements)
        {
            if(!genomeRegion.chromosome().equals(svData.chromosome(useStart)))
                continue;

            // test if the SV falls within the LE +/- a buffer
            if(svData.position(useStart) >= genomeRegion.start() - PERMITTED_DISTANCE
                    && svData.position(useStart) <= genomeRegion.end() + PERMITTED_DISTANCE)
            {
                LOGGER.debug("var({}) found in identified line element({}.>{})",
                        svData.posId(), genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end());
                return IDENTIFIED_LINE_ELEMENT;
            }
        }

        return NO_LINE_ELEMENT;
    }

}
