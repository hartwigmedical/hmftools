package com.hartwig.hmftools.telo.analysers;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.telo.ReadGroup;
import com.hartwig.hmftools.telo.TeloUtils;

// analyse the telo bam and come up with our own score
// of F1, F2, F4
public class TelomereReadsAnalyser
{
    private final Map<ReadGroup.FragmentType, Integer> mFragmentTypeCount = new HashMap<>();

    public int getFragmentTypeCount(ReadGroup.FragmentType fragType)
    {
        return mFragmentTypeCount.getOrDefault(fragType, 0);
    }

    public void onReadGroup(final ReadGroup readGroup)
    {
        // we want to classify this read group
        ReadGroup.FragmentType fragType = TeloUtils.classifyFragment(readGroup);
        mFragmentTypeCount.put(fragType, mFragmentTypeCount.getOrDefault(fragType, 0) + 1);
    }
}
