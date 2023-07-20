package com.hartwig.hmftools.ctdna.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class SampleData
{
    public final String PatientId;
    public final String TumorId;
    public final List<String> CtDnaSamples;
    public final String VcfTag;

    public SampleData(final String patientId, final String tumorId, final List<String> ctDnaSamples, final String vcfTag)
    {
        PatientId = patientId;
        TumorId = tumorId;
        CtDnaSamples = ctDnaSamples;
        VcfTag = vcfTag;
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(", ");
        CtDnaSamples.forEach(x -> sj.add(x));
        return format("patient(%s) tumor(%s) ctDnaSamples(%s)", PatientId, TumorId, sj);
    }

    public static List<String> ctDnaSamplesFromStr(final String ctDnaSamples)
    {
        return Arrays.stream(ctDnaSamples.split(ITEM_DELIM, -1)).collect(Collectors.toList());
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
            int ctdnaIndex = fieldsIndexMap.get("CtDnaSampleIds");
            Integer vcfIndex = fieldsIndexMap.get("VcfTag");

            for(String line : fileContents)
            {
                if(line.startsWith("#") || line.isEmpty())
                    continue;

                String[] values = line.split(CSV_DELIM, -1);
                String vcfTag = vcfIndex != null && vcfIndex < values.length ? values[vcfIndex] : "";

                samples.add(new SampleData(
                        values[patientIndex], values[tumorIndex], ctDnaSamplesFromStr(values[ctdnaIndex]), vcfTag));
            }
        }
        catch (IOException e)
        {
            CT_LOGGER.error("failed to read sample data file({}): {}", filename, e.toString());
            return null;
        }

        return samples;
    }
}
