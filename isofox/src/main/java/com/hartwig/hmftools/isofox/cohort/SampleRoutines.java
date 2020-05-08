package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class SampleRoutines
{
    private final CohortConfig mConfig;

    public SampleRoutines(final CohortConfig config)
    {
        mConfig = config;
    }

    public void convertMutationCohortFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("no file({}) found", filename);
            return;
        }

        try
        {
            final String outputFileName = mConfig.formCohortFilename("sample_mut_vs_wt_data.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("Cohort,CancerType,SampleId,GeneName");
            writer.newLine();

            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final String[] fieldNames = lines.get(0).split("\t");
            lines.remove(0);

            // Entity	Gene	Group1_Mut	Group1_Wt	Group2_Mut	Group2_Wt

            for(final String data : lines)
            {
                final String[] items = data.split("\t");

                final String cancerType = items[0];
                final String geneName = items[1];

                for(int i = 0; i < 4; ++i)
                {
                    final String mutData = items[i + 2].replaceAll("\"", "");

                    if(mutData.isEmpty())
                        continue;

                    final String[] samples = mutData.split(",", -1);
                    final String cohortName = fieldNames[i + 2];

                    for(int j = 0; j < samples.length; ++j)
                    {
                        writer.write(String.format("%s,%s,%s,%s", cohortName, cancerType, samples[j], geneName));
                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to convert sample mutations file({}): {}", filename, e.toString());
        }
    }
}
