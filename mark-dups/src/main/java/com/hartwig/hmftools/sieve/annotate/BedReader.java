package com.hartwig.hmftools.sieve.annotate;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public class BedReader
{
    @NotNull
    static public List<AnnotatedBedRecord> readFromFile(@NotNull final String filepath)
    {
        final List<AnnotatedBedRecord> records = new ArrayList<>();

        try
        {
            final BufferedReader reader = new BufferedReader(new FileReader(filepath));

            // Drop header.
            reader.readLine();

            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] fields = line.split("\t");
                if(fields.length != 6)
                {
                    MD_LOGGER.error("BED file {} contains a record with {} fields. There should only be six fields. ", filepath, fields.length);
                    System.exit(1);
                }

                final AnnotatedBedRecord record = makeBedRecord(fields);
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

    @NotNull
    static private AnnotatedBedRecord makeBedRecord(@NotNull final String[] fields)
    {
        final String chromosome = fields[0];

        int posStart = 0;
        try
        {
            posStart = Integer.parseInt(fields[1]);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("Failed to parse 'posStart' field to an int while reading a BED file.");
            System.exit(1);
        }

        int posEnd = 0;
        try
        {
            posEnd = Integer.parseInt(fields[2]);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("Failed to parse 'posEnd' field to an int while reading a BED file.");
            System.exit(1);
        }

        long sampleCount = 0;
        try
        {
            sampleCount = Long.parseLong(fields[3]);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("Failed to parse 'sampleCount' field to a long while reading a BED file.");
            System.exit(1);
        }

        long depthMin = 0;
        try
        {
            depthMin = Long.parseLong(fields[4]);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("Failed to parse 'depthMin' field to a long while reading a BED file.");
            System.exit(1);
        }

        long depthMax = 0;
        try
        {
            depthMax = Long.parseLong(fields[5]);
        }
        catch(NumberFormatException e)
        {
            MD_LOGGER.error("Failed to parse 'depthMax' field to a long while reading a BED file.");
            System.exit(1);
        }

        return new AnnotatedBedRecord(chromosome, posStart, posEnd, sampleCount, depthMin, depthMax);
    }
}
