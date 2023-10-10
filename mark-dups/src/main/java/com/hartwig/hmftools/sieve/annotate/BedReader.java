package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

// TODO(m_cooper): Code duplication?
public class BedReader
{
    private static final String DELIMITER = "\t";
    private static final int NUM_FIELDS = 6;

    static public List<BlacklistRegion> readFromFile(final String filepath)
    {
        final List<BlacklistRegion> records = new ArrayList<>();

        try
        {
            final BufferedReader reader = new BufferedReader(new FileReader(filepath));

            // Drop header.
            reader.readLine();

            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] fields = line.split(DELIMITER);
                if(fields.length != NUM_FIELDS)
                {
                    MD_LOGGER.error("BED file {} contains a record with {} fields. There should only be {} fields. ", filepath, fields.length, NUM_FIELDS);
                    System.exit(1);
                }

                final BlacklistRegion record = parseFields(fields, filepath);
                records.add(record);
            }

            reader.close();

        }
        catch(Exception e)
        {
            MD_LOGGER.error("An exception was raised while reading the BED file {}: {}", filepath, e.toString());
            System.exit(1);
        }

        return records;
    }

    static private BlacklistRegion parseFields(final String[] fields, final String filepath)
    {
        String chromosome = fields[0];
        int posStart = parseInt(fields[1], "PosStart", filepath);
        int posEnd = parseInt(fields[2], "PosEnd", filepath);
        int sampleCount = parseInt(fields[3], "SampleCount", filepath);
        int depthMin = parseInt(fields[4], "DepthMin", filepath);
        int depthMax = parseInt(fields[5], "DepthMax", filepath);

        return new BlacklistRegion(chromosome, posStart, posEnd, sampleCount, depthMin, depthMax);
    }

    static private int parseInt(final String str, final String fieldName, final String filepath)
    {
        try
        {
            return Integer.parseInt(str);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("While reading a record in the BED file {}, failed to parse the '{}' field to an int.", filepath, fieldName);
            System.exit(1);
        }

        return 0;
    }
}
