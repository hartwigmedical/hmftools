package com.hartwig.hmftools.neo.scorer;

import static com.hartwig.hmftools.common.utils.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.scorer.NeoScorerConfig.SAMPLE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;

public class SampleData
{
    public final String Id;
    public final String CancerType;
    public final boolean HasRna;

    public SampleData(final String id, final String cancerType, boolean hasRna)
    {
        Id = id;
        CancerType = cancerType;
        HasRna = hasRna;
    }

    public static List<SampleData> loadFromConfig(final CommandLine cmd)
    {
        List<SampleData> samples = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            String filename = cmd.getOptionValue(SAMPLE_ID_FILE);
            try
            {
                final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

                if(fileContents.isEmpty())
                    return samples;

                String header = fileContents.get(0);

                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
                fileContents.remove(0);

                int sampleIdIndex = fieldsIndexMap.get("SampleId");
                Integer cancerTypeIndex = fieldsIndexMap.get("CancerType");
                Integer rnaIndex = fieldsIndexMap.get("HasRna");

                for(String line : fileContents)
                {
                    if(line.startsWith("#") || line.isEmpty())
                        continue;

                    String[] values = line.split(CSV_DELIM, -1);
                    String sampleId = values[sampleIdIndex];
                    String cancerType = cancerTypeIndex != null ? values[cancerTypeIndex] : "";
                    boolean hasRna = rnaIndex != null ? Boolean.parseBoolean(values[rnaIndex]) : false;
                    samples.add(new SampleData(sampleId, cancerType, hasRna));
                }
            }
            catch (IOException e)
            {
                NE_LOGGER.error("failed to read sample ID file({}): {}", filename, e.toString());
            }
        }
        else if(cmd.hasOption(SAMPLE))
        {
            String sampleData = cmd.getOptionValue(SAMPLE);
            String[] values = sampleData.split(CSV_DELIM, -1);
            String sampleId = values[0];
            String cancerType = values.length >= 2 ?  values[1] : "";
            boolean hasRna = values.length >= 3 ? Boolean.parseBoolean(values[2]) : false;
            samples.add(new SampleData(sampleId, cancerType, hasRna));
        }

        return samples;
    }
}
