package com.hartwig.hmftools.neo;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.neo.epitope.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class NeoCommon
{
    public static final int DOWNSTREAM_PRE_GENE_DISTANCE = 100000; // in concordance with Linx

    public static final String NEO_FILE_ID = ".neo.";

    public static final String OUTPUT_ID = "output_id";
    public static final String THREADS = "threads";

    public static final Logger NE_LOGGER = LogManager.getLogger(NeoCommon.class);

    public static void loadSampleIdsFile(final String filename, final List<String> sampleIds)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.warn("invalid sampleId file({})", filename);
            return;
        }

        sampleIds.addAll(ConfigUtils.loadSampleIdsFile(filename));
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
