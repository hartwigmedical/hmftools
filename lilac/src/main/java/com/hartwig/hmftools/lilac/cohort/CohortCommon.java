package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

public class CohortCommon
{
    public static final String SAMPLE_IDS_FILE = "sample_ids_file";
    public static final String SAMPLE_FILES_DIR = "sample_files_dir";

    public static List<String> loadSampleIds(final String fileName)
    {
        List<String> sampleIds = Lists.newArrayList();

        try
        {
            final List<String> lines = Files.readAllLines(new File(fileName).toPath());
            lines.stream().filter(x -> !x.equals("SampleId")).forEach(x -> sampleIds.add(x));
            LL_LOGGER.info("load {} samples", sampleIds.size());
        }
        catch (IOException e)
        {
            LL_LOGGER.warn("failed to load sample file({}): {}", fileName, e.toString());
        }

        return sampleIds;

    }
}
