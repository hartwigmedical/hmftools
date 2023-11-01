package com.hartwig.hmftools.cobalt.targeted;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class TargetRegionEnrichment
{
    private final List<GenomePosition> mTargetedRegions = new ArrayList<>();
    private final Map<GenomePosition, Double> mTargetRelativeEnrichment = new TreeMap<>(GenomePosition::compare);

    public List<GenomePosition> getTargetedRegions()
    {
        return mTargetedRegions;
    }

    public Map<GenomePosition, Double> getTargetRelativeEnrichment()
    {
        return mTargetRelativeEnrichment;
    }

    private static final char DELIMITER = '\t';

    public static TargetRegionEnrichment fromTsv(String fileName) throws IOException
    {
        TargetRegionEnrichment targetRegionEnrichment = new TargetRegionEnrichment();

        try (BufferedReader reader = new BufferedReader(new FileReader(fileName)))
        {
            CSVFormat format = CSVFormat.Builder.create()
                    .setDelimiter(DELIMITER)
                    .setRecordSeparator('\n')
                    .setHeader().setSkipHeaderRecord(true)
                    .build();
            Iterable<CSVRecord> records = format.parse(reader);

            for (CSVRecord record : records)
            {
                String chromosome = record.get("chromosome").intern();
                int position = (int) Double.parseDouble(record.get("position"));
                Double relativeEnrichment = parseDoubleOrNull(record.get("relativeEnrichment"));
                GenomePosition genomePosition = GenomePositions.create(chromosome, position);
                targetRegionEnrichment.mTargetedRegions.add(genomePosition);
                if (relativeEnrichment != null && !Double.isNaN(relativeEnrichment))
                {
                    targetRegionEnrichment.mTargetRelativeEnrichment.put(genomePosition, relativeEnrichment);
                }
            }
        }

        return targetRegionEnrichment;
    }

    private static Double parseDoubleOrNull(String value)
    {
        try
        {
            return Double.parseDouble(value);
        }
        catch (NumberFormatException e)
        {
            return null;
        }
    }
}
