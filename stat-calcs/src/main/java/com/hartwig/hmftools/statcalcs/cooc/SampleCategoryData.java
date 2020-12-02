package com.hartwig.hmftools.statcalcs.cooc;

import java.util.List;

import com.google.common.collect.Lists;

public class SampleCategoryData
{
    public static final int SAMPLE_CAT_1_INDEX = 0;
    public static final int SAMPLE_CAT_2_INDEX = 1;

    public final String SampleId;

    private final List<String[]> mCategoryData;

    public SampleCategoryData(final String sampleId)
    {
        SampleId = sampleId;
        mCategoryData = Lists.newArrayList();
    }

    public void addCategoryData(final String cat1, final String cat2)
    {
        String[] values = {cat1, cat2};
        mCategoryData.add(values);
    }

    public final List<String[]> getCategoryData() { return mCategoryData; }
}
