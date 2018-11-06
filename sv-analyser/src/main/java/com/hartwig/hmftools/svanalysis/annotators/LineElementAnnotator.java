package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LineElementAnnotator {

    public static String KNOWN_LINE_ELEMENT = "Known";
    public static String IDENTIFIED_LINE_ELEMENT = "Ident";
    public static String NO_LINE_ELEMENT = "None";
    public static String SUSPECTED_LINE_ELEMENT = "Suspect";

    private static String CSV_LE_TYPE_IDENTIFIED = "Identified";

    private List<GenomeRegion> mKnownLineElements;
    private List<GenomeRegion> mIdentifiedLineElements;
    private static int PERMITTED_DISTANCE = 5000;

    private static int POLY_A_T_LENGTH = 5;

    private String mPolyAMotif;
    private String mPolyTMotif;

    private static final Logger LOGGER = LogManager.getLogger(FragileSiteAnnotator.class);

    public LineElementAnnotator()
    {
        mKnownLineElements = Lists.newArrayList();
        mIdentifiedLineElements = Lists.newArrayList();

        mPolyAMotif = "";
        mPolyTMotif = "";
        for(int i = 0; i < POLY_A_T_LENGTH; ++i)
        {
            mPolyTMotif += "T";
            mPolyAMotif += "A";
        }
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

    public String isLineElement(final SvVarData svData, final boolean useStart)
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

    public void setSuspectedLineElements(final Map<String, List<SvBreakend>> chrBreakendMap, int proximityLength)
    {
        // if there are 3 or more BNDs within the standard proximity window and they connect to 2 or more other arms,
        // classify these as suspected LINE elements
        for(Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            // final String chromosome = entry.getKey();

            List<SvBreakend> breakendList = entry.getValue();

            List<Long> positions = Lists.newArrayList();
            List<SvBreakend> potentialLineSVs = Lists.newArrayList();

            boolean isLineGroup = false;

            for (final SvBreakend breakend : breakendList)
            {
                long newPosition = breakend.position();

                // remove earlier breakends which are now too far from the current one
                int index = 0;
                while(index < positions.size())
                {
                    if(newPosition - positions.get(0) <= proximityLength)
                        break;

                    positions.remove(index);
                    potentialLineSVs.remove(index);
                    isLineGroup = false; // require reassessment below
                }

                positions.add(newPosition);
                potentialLineSVs.add(breakend);

                if(positions.size() < 3)
                {
                    isLineGroup = false;
                    continue;
                }

                if(!isLineGroup)
                {
                    boolean hasMultipleRemoteArms = false;
                    int bndCount = 0;
                    boolean hasPolyATMotify = false;

                    String currentOtherChr = "";
                    for (final SvBreakend lineBreakend : potentialLineSVs)
                    {
                        final SvVarData var = breakend.getSV();

                        if (var.type() == BND)
                        {
                            ++bndCount;

                            if(var.getSvData().insertSequence().contains(mPolyAMotif) || var.getSvData().insertSequence().contains(mPolyTMotif))
                                hasPolyATMotify = true;

                            final String otherChr = var.chromosome(!lineBreakend.usesStart());

                            if (currentOtherChr.isEmpty())
                            {
                                currentOtherChr = otherChr;
                            }
                            else if (currentOtherChr.equals(otherChr))
                            {
                                continue;
                            }
                            else
                            {
                                hasMultipleRemoteArms = true;
                            }
                        }

                        if (bndCount >= 3 && hasMultipleRemoteArms && hasPolyATMotify)
                        {
                            isLineGroup = true;
                            break;
                        }
                    }
                }

                if(isLineGroup)
                {
                    for (SvBreakend lineBreakend : potentialLineSVs)
                    {
                        if(!lineBreakend.getSV().isLineElement(lineBreakend.usesStart()))
                            lineBreakend.getSV().setLineElement(SUSPECTED_LINE_ELEMENT, lineBreakend.usesStart());
                    }
                }
            }
        }

    }

}
