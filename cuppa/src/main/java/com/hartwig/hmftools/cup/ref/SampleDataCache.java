package com.hartwig.hmftools.cup.ref;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.utils.CuppaConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.utils.CuppaConstants.DATA_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

@Deprecated
public class SampleDataCache
{
    public final List<String> SampleIds;
    public final List<SampleData> SampleDataList;
    public SampleData SpecificSample;

    public final Map<String,List<SampleData>> RefCancerSampleData; // map of known cancer types to their ref samples
    public final Map<String,String> RefSampleCancerTypeMap; // map of ref sample to cancer type
    public final List<SampleData> RefSampleDataList; // includes type 'Other' if in ref data
    public final Map<String,Integer> RefSampleRnaReadLength;

    public final Map<String,String> RefCancerMappings; // subtypes to parent types

    private boolean mIsValid;

    public SampleDataCache()
    {
        SampleIds = Lists.newArrayList();
        SampleDataList = Lists.newArrayList();

        RefCancerSampleData = Maps.newHashMap();
        RefSampleCancerTypeMap = Maps.newHashMap();
        RefSampleDataList = Lists.newArrayList();
        RefSampleRnaReadLength = Maps.newHashMap();
        RefCancerMappings = Maps.newHashMap();

        SpecificSample = null;
        mIsValid = true;
    }

    public boolean isValid() { return mIsValid; }

    public void addRefSample(final SampleData sample)
    {
        if(sample.isKnownCancerType())
        {
            List<SampleData> cancerSampleData = RefCancerSampleData.get(sample.cancerType());

            if(cancerSampleData == null)
                RefCancerSampleData.put(sample.cancerType(), Lists.newArrayList(sample));
            else
                cancerSampleData.add(sample);

            if(!sample.cancerType().equals(sample.cancerMainType()))
            {
                RefCancerMappings.put(sample.cancerType(), sample.cancerMainType());
            }
        }

        RefSampleCancerTypeMap.put(sample.Id, sample.cancerType());
        RefSampleDataList.add(sample);

        if(sample.rnaReadLength() > 0)
            RefSampleRnaReadLength.put(sample.Id, sample.rnaReadLength());
    }

    public void loadReferenceSampleData(final String refSampleDataFile)
    {
        if(refSampleDataFile == null)
        {
            mIsValid = false;
            return;
        }

        try
        {
            final List<String> fileData = Files.readAllLines(new File(refSampleDataFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SampleData sample = SampleData.from(fieldsIndexMap, line);
                addRefSample(sample);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read reference sample data file({}): {}", refSampleDataFile, e.toString());
            mIsValid = false;
        }
    }
}
