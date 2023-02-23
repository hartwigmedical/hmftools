package com.hartwig.hmftools.fastqtools;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;

import java.util.Map;

import com.google.common.collect.Maps;

public class RepeatFrequencies
{
    private final int[] mBaseCounts;
    private final Map<String,Integer> mRepeatFrequencies;

    public RepeatFrequencies()
    {
        mRepeatFrequencies = Maps.newHashMap();
        mBaseCounts = new int[DNA_BASES.length];
    }

    public int[] baseCounts() { return mBaseCounts; }
    public Map<String,Integer> repeatFrequencies() { return mRepeatFrequencies; }

    public void registerBase(final char base)
    {
        for(int i = 0; i < mBaseCounts.length; ++i)
        {
            if(DNA_BASES[i] == base)
            {
                ++mBaseCounts[i];
                return;
            }
        }
    }

    public void registerSequence(final String repeatSequence)
    {
        if(repeatSequence.length() < 2)
            return;

        Integer count = mRepeatFrequencies.get(repeatSequence);
        mRepeatFrequencies.put(repeatSequence, count != null ? count + 1 : 1);
    }
}
