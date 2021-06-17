package com.hartwig.hmftools.neo;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NeoCommon
{
    public static final int DOWNSTREAM_PRE_GENE_DISTANCE = 100000; // in concordance with Linx

    public static final String IM_FILE_ID = ".imu.";
    public static final String LOG_DEBUG = "log_debug";

    public static final Logger NE_LOGGER = LogManager.getLogger(NeoCommon.class);

    public static void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("GeneId"))
            {
                // check for header row
                fileContents.remove(0);
            }

            geneIdList.addAll(fileContents.stream()
                    .filter(x -> !x.isEmpty())
                    .filter(x -> !x.contains("GeneId"))
                    .filter(x -> !x.startsWith("#"))
                    .map(x -> x.split(",")[0]).collect(Collectors.toList()));
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

    public static void loadSampleIdsFile(final String filename, final List<String> sampleIds)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("invalid sampleId file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("SampleId"))
                fileContents.remove(0);

            for(String sampleData : fileContents)
            {
                if(sampleData.contains(DELIMITER))
                {
                    String[] items = sampleData.split(DELIMITER, -1);
                    sampleIds.add(items[0]);
                }
                else
                {
                    sampleIds.add(sampleData);
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load sampleId file({}): {}", filename, e.toString());
        }
    }

    public static void loadSampleDataFile(final String filename, final List<SampleData> samples)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("invalid sampleId file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("SampleId"))
                fileContents.remove(0);

            for(String data : fileContents)
            {
                SampleData sample = SampleData.fromCsv(data);

                if(sample == null)
                {
                    NE_LOGGER.error("invalid sample data record: {}", data);
                    continue;
                }

                samples.add(sample);
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load sample data file({}): {}", filename, e.toString());
        }
    }

}
