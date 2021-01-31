package com.hartwig.hmftools.imuno.common;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.imuno.neo.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ImunoCommon
{
    public static final int DOWNSTREAM_PRE_GENE_DISTANCE = 100000; // in concordance with Linx

    public static final String LOG_DEBUG = "log_debug";

    public static final Logger IM_LOGGER = LogManager.getLogger(ImunoCommon.class);

    public static void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            IM_LOGGER.warn("invalid gene ID file({})", filename);
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
            IM_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

    public static void loadSampleIdsFile(final String filename, final List<String> sampleIds)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            IM_LOGGER.warn("invalid sampleId file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("SampleId"))
                fileContents.remove(0);

            sampleIds.addAll(fileContents);
        }
        catch (IOException e)
        {
            IM_LOGGER.warn("failed to load sampleId file({}): {}", filename, e.toString());
        }
    }

    public static void loadSampleDataFile(final String filename, final List<SampleData> samples)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            IM_LOGGER.warn("invalid sampleId file({})", filename);
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
                    IM_LOGGER.error("invalid sample data record: {}", data);
                    continue;
                }

                samples.add(sample);
            }
        }
        catch (IOException e)
        {
            IM_LOGGER.warn("failed to load sample data file({}): {}", filename, e.toString());
        }
    }

}
