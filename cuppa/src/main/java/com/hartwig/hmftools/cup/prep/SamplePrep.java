package com.hartwig.hmftools.cup.prep;

import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;

public class SamplePrep implements Callable
{
    private final List<String> mSampleIds;

    public SamplePrep()
    {
        mSampleIds = Lists.newArrayList();
    }

    public void addSampleId(final String sampleId) { mSampleIds.add(sampleId); }

    @Override
    public Long call() throws Exception
    {
        for(String sampleId : mSampleIds)
        {

        }

        return (long)0;
    }

    public void prepareSampleData(final String sampleId)
    {

    }

}
