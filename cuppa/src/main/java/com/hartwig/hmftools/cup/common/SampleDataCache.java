package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.CuppaConfig.CANCER_SUBTYPE_OTHER;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_UNKNOWN;
import static com.hartwig.hmftools.cup.common.SampleData.RNA_READ_LENGTH_NONE;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.traits.SampleTraitsData;

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
    public final Map<String,SampleTraitsData> SampleTraitsData;

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
        SampleTraitsData = Maps.newHashMap();
        RefCancerMappings = Maps.newHashMap();

        SpecificSample = null;
        mIsValid = true;
    }

    public boolean isValid() { return mIsValid; }
    public boolean isSingleSample() { return SpecificSample != null; }
    public boolean isMultiSample() { return !isSingleSample(); }
    public boolean isMultiSampleNonRef() { return SampleIds.size() > 1 && RefSampleCancerTypeMap.size() != SampleIds.size(); }

    public int getCancerSampleCount(final String cancerType)
    {
        final List<SampleData> samples = RefCancerSampleData.get(cancerType);
        return samples != null ? samples.size() : 0;
    }

    public boolean hasRefSample(final String sampleId) { return RefSampleCancerTypeMap.containsKey(sampleId); }
    public boolean hasRefCancerType(final String cancerType) { return RefCancerSampleData.containsKey(cancerType); }

    public final List<String> refSampleIds(boolean onlyKnownTypes)
    {
        return RefSampleDataList.stream()
                .filter(x -> !onlyKnownTypes || x.isKnownCancerType())
                .map(x -> x.Id).collect(Collectors.toList());
    }

    public SampleData findRefSampleData(final String sampleId)
    {
        for(List<SampleData> sampleDataList : RefCancerSampleData.values())
        {
            final SampleData refSample = sampleDataList.stream().filter(x -> x.Id.equals(sampleId)).findFirst().orElse(null);
            if(refSample != null)
                return refSample;
        }

        return null;
    }

    public SampleData findSampleData(final String sampleId)
    {
        return SampleDataList.stream().filter(x -> x.Id.equals(sampleId)).findFirst().orElse(null);
    }

    public int getRefSampleRnaReadLength(final String sampleId)
    {
        Integer readLength = RefSampleRnaReadLength.get(sampleId);
        return readLength != null ? readLength : RNA_READ_LENGTH_NONE;
    }

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

    public SampleData addTestSample(final SampleData sample)
    {
        if(SampleIds.contains(sample.Id))
        {
            CUP_LOGGER.error("attempted to add sample({}) twice", sample.Id);
            return null;
        }

        SampleData testSample = sample;
        String refCancerType = RefSampleCancerTypeMap.get(sample.Id);

        if(refCancerType != null)
        {
            // override if known
            testSample = findRefSampleData(sample.Id);

            if(testSample == null)
            {
                CUP_LOGGER.error("sample({} ct={}) marked as ref but not found", sample.Id, sample.cancerType());
                return null;
            }
            else
            {
                testSample.setRefSample();
                testSample.setCancerType(refCancerType);
            }
        }

        SampleDataList.add(testSample);
        SampleIds.add(testSample.Id);
        return testSample;
    }

    public void loadSampleData(final String specificSampleId, final int sampleRnaReadLength, final String sampleDataFile)
    {
        if(specificSampleId != null)
        {
            SampleData sample = new SampleData(specificSampleId, "", "");
            sample.setRnaReadLength(sampleRnaReadLength);
            SpecificSample = addTestSample(sample);

            if(SpecificSample.isRefSample())
            {
                CUP_LOGGER.info("sample({}) cancerType({}) is a ref sample", SpecificSample.Id, SpecificSample.cancerType());
            }
        }
        else if(sampleDataFile != null)
        {
            try
            {
                final List<String> fileData = Files.readAllLines(new File(sampleDataFile).toPath());

                final String header = fileData.get(0);
                fileData.remove(0);

                final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

                for(final String line : fileData)
                {
                    SampleData sample = SampleData.from(fieldsIndexMap, line);
                    addTestSample(sample);
                }
            }
            catch (IOException e)
            {
                CUP_LOGGER.error("failed to read sample data file({}): {}", sampleDataFile, e.toString());
                mIsValid = false;
            }
        }
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
