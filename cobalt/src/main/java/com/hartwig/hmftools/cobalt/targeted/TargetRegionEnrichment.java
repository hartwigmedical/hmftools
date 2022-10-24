package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class TargetRegionEnrichment
{
    private final List<GenomePosition> mTargetedRegions;
    private final Map<GenomePosition, Double> mTargetRelativeEnrichment;

    private static final String DELIM = "\t";

    public TargetRegionEnrichment(final List<GenomePosition> targetedRegions, final Map<GenomePosition, Double> targetRelativeEnrichment)
    {
        mTargetedRegions = targetedRegions;
        mTargetRelativeEnrichment = targetRelativeEnrichment;
        // mTargetRelativeEnrichment: MutableMap<GenomePosition, Double> = TreeMap(GenomePosition::compare);
    }

    public List<GenomePosition> regions() { return mTargetedRegions; }
    public Map<GenomePosition,Double> regionEnrichment() { return mTargetRelativeEnrichment; }

    public static TargetRegionEnrichment load(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);

            int chrIndex = fieldsIndexMap.get("chromosome");
            int posIndex = fieldsIndexMap.get("position");
            int reIndex = fieldsIndexMap.get("relativeEnrichment");

            lines.remove(0);

            List<GenomePosition> targetedRegions = Lists.newArrayList();
            Map<GenomePosition,Double> targetRelativeEnrichment = Maps.newTreeMap(GenomePosition::compare);
            // Map<GenomePosition, Double> targetRelativeEnrichment = Maps.newHashMap();

            for(String line : lines)
            {
                String[] values = line.split(DELIM, -1);

                String chromosome = values[chrIndex];
                int position = Integer.parseInt(values[posIndex]);

                GenomePosition genomePosition = GenomePositions.create(chromosome, position);

                targetedRegions.add(genomePosition);

                try
                {
                    double relativeEnrichment = Double.parseDouble(values[reIndex]);

                    if(!Double.isNaN(relativeEnrichment))
                        targetRelativeEnrichment.put(genomePosition, relativeEnrichment);
                }
                catch(Exception e)
                {}
            }

            return new TargetRegionEnrichment(targetedRegions, targetRelativeEnrichment);
        }
        catch(Exception e)
        {
            CB_LOGGER.error("failed to load target enrichment regions: {}", e.toString());
            return null;
        }
    }
}