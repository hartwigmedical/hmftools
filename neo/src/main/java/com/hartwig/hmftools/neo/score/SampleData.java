package com.hartwig.hmftools.neo.score;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.IGNORE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CANCER_TYPE;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.score.NeoScorerConfig.CANCER_TYPE;
import static com.hartwig.hmftools.neo.score.NeoScorerConfig.RNA_SOMATIC_VCF;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class SampleData
{
    public final String TumorId;
    public final String RnaSampleId;
    public final String CancerType;
    public final boolean HasRna;

    public SampleData(final String tumorId, final String rnaSampleId, final String cancerType, boolean hasRna)
    {
        TumorId = tumorId;
        RnaSampleId = rnaSampleId;
        CancerType = cancerType;
        HasRna = hasRna;
    }

    public static List<SampleData> loadFromConfig(final ConfigBuilder configBuilder)
    {
        List<SampleData> samples = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            String filename = configBuilder.getValue(SAMPLE_ID_FILE);
            try
            {
                final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

                if(fileContents.isEmpty())
                    return samples;

                String header = fileContents.get(0);

                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
                fileContents.remove(0);

                int sampleIdIndex = fieldsIndexMap.get(FLD_SAMPLE_ID);
                Integer cancerTypeIndex = fieldsIndexMap.get(FLD_CANCER_TYPE);
                Integer rnaIndex = fieldsIndexMap.get("HasRna");

                for(String line : fileContents)
                {
                    if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                        continue;

                    String[] values = line.split(CSV_DELIM, -1);
                    String sampleId = values[sampleIdIndex];
                    String cancerType = cancerTypeIndex != null ? values[cancerTypeIndex] : "";
                    boolean hasRna = rnaIndex != null ? Boolean.parseBoolean(values[rnaIndex]) : false;
                    samples.add(new SampleData(sampleId, null, cancerType, hasRna));
                }
            }
            catch (IOException e)
            {
                NE_LOGGER.error("failed to read sample ID file({}): {}", filename, e.toString());
            }
        }
        else if(configBuilder.hasValue(SAMPLE))
        {
            String sampleId = configBuilder.getValue(SAMPLE);
            String rnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
            String cancerType = configBuilder.getValue(CANCER_TYPE, "");
            boolean hasRna = configBuilder.hasValue(ISOFOX_DIR_CFG) && configBuilder.hasValue(RNA_SOMATIC_VCF);
            samples.add(new SampleData(sampleId, rnaSampleId, cancerType, hasRna));
        }

        return samples;
    }
}
