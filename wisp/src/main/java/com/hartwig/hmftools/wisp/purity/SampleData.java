package com.hartwig.hmftools.wisp.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.IGNORE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.wisp.common.CommonUtils.BATCH_CONTROL_TAG;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.wisp.common.CommonUtils;

public class SampleData
{
    public final String PatientId;
    public final String TumorId;
    public final String AmberExtraTumorId;
    public final List<String> SampleIds;
    public final String VcfTag;
    public final boolean IsPanel;

    public SampleData(
            final String patientId, final String tumorId, final List<String> sampleIds, final String vcfTag, final boolean isPanel,
            final String amberExtraTumorId)
    {
        PatientId = patientId;
        TumorId = tumorId;
        SampleIds = sampleIds;
        VcfTag = vcfTag;
        IsPanel = isPanel;
        AmberExtraTumorId = amberExtraTumorId;
    }

    public boolean isBatchControl() { return VcfTag != null && VcfTag.contains(BATCH_CONTROL_TAG); }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(", ");
        SampleIds.forEach(x -> sj.add(x));
        return format("patient(%s) tumor(%s) samples(%s)", PatientId, TumorId, sj);
    }

    public static List<String> sampleIdsFromStr(final String sampleIds)
    {
        return Arrays.stream(sampleIds.split(ITEM_DELIM, -1)).collect(Collectors.toList());
    }

    public static List<SampleData> loadSampleDataFile(final String filename)
    {
        List<SampleData> samples = Lists.newArrayList();

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            String header = fileContents.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
            fileContents.remove(0);

            int patientIndex = fieldsIndexMap.get("PatientId");
            int tumorIndex = fieldsIndexMap.get("TumorId");
            int sampleIdsIndex = fieldsIndexMap.get("SampleIds");
            Integer vcfIndex = fieldsIndexMap.get("VcfTag");
            Integer isPanelIndex = fieldsIndexMap.get("IsPanel");
            Integer amberExtraTumorIdIndex = fieldsIndexMap.get("AmberExtraTumorId");

            for(String line : fileContents)
            {
                if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                    continue;

                String[] values = line.split(CSV_DELIM, -1);
                String vcfTag = vcfIndex != null && vcfIndex < values.length ? values[vcfIndex] : "";

                boolean isPanel = isPanelIndex != null && isPanelIndex < values.length ?
                        Boolean.parseBoolean(values[isPanelIndex]) : false;

                String patientId = values[patientIndex];
                String tumorId = values[tumorIndex];
                String amberExtraTumorId = amberExtraTumorIdIndex != null ? values[amberExtraTumorIdIndex] : null;
                List<String> sampleIds = sampleIdsFromStr(values[sampleIdsIndex]);

                samples.add(new SampleData(patientId, tumorId, sampleIds, vcfTag, isPanel, amberExtraTumorId));
            }
        }
        catch (IOException e)
        {
            CommonUtils.CT_LOGGER.error("failed to read sample data file({}): {}", filename, e.toString());
            return null;
        }

        return samples;
    }
}
