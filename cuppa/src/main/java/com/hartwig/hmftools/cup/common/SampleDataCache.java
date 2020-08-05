package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SAMPLE_DATA_FILE;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;

public class SampleDataCache
{
    public final List<SampleData> SampleDataList;
    public SampleData SpecificSample;
    public final Map<String,List<SampleData>> CancerSampleData;
    public final Map<String,String> SampleCancerTypeMap;

    public SampleDataCache()
    {
        SampleDataList = Lists.newArrayList();
        CancerSampleData = Maps.newHashMap();
        SampleCancerTypeMap = Maps.newHashMap();
        SpecificSample = null;
    }

    public void loadSampleData(final String specificSampleData, final String sampleDataFile)
    {
        if(specificSampleData != null)
        {
            final String[] sampleItems = specificSampleData.split(SUBSET_DELIM, -1);
            String sampleId = sampleItems[0];
            String cancerType = CANCER_TYPE_UNKNOWN;
            String cancerSubtype = CANCER_SUBTYPE_OTHER;

            if(sampleItems.length >= 2)
                cancerType = sampleItems[1];

            if(sampleItems.length == 3)
                cancerSubtype = sampleItems[2];

            SpecificSample = new SampleData(sampleId, cancerType, cancerSubtype);
            SampleDataList.add(SpecificSample);
        }

        if(sampleDataFile != null)
        {
            try
            {
                final List<String> fileData = Files.readAllLines(new File(sampleDataFile).toPath());

                final String header = fileData.get(0);
                fileData.remove(0);

                final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

                for(final String line : fileData)
                {
                    final String[] items = line.split(DATA_DELIM, -1);
                    final String sampleId = items[fieldsIndexMap.get("SampleId")];
                    final String cancerType = items[fieldsIndexMap.get("CancerType")];

                    final String cancerSubtype = fieldsIndexMap.containsKey("CancerSubtype") ?
                            items[fieldsIndexMap.get("CancerSubtype")] : CANCER_SUBTYPE_OTHER;

                    if(SpecificSample != null && SpecificSample.Id.equals(sampleId))
                        continue;

                    SampleDataList.add(new SampleData(sampleId, cancerType, cancerSubtype));
                }
            }
            catch (IOException e)
            {
                CUP_LOGGER.error("failed to read sample data file({}): {}", sampleDataFile, e.toString());
            }
        }

        // build a cache of samples per cancer type
        for(SampleData sampleData : SampleDataList)
        {
            SampleCancerTypeMap.put(sampleData.Id, sampleData.CancerType);

            List<SampleData> cancerSampleData = CancerSampleData.get(sampleData.CancerType);

            if(cancerSampleData == null)
                CancerSampleData.put(sampleData.CancerType, Lists.newArrayList(sampleData));
            else
                cancerSampleData.add(sampleData);
        }
    }
}
