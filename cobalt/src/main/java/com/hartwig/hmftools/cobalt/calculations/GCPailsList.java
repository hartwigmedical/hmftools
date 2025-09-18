package com.hartwig.hmftools.cobalt.calculations;

import java.util.ArrayList;
import java.util.List;

public class GCPailsList
{
    private final List<GCPail> mBuckets = new ArrayList<>(101);

    public GCPailsList()
    {
        for (int i = 0; i < 101; i++)
        {
            mBuckets.add(new GCPail(i));
        }
    }

    public GCPail getGCPail(double gcContent)
    {
        return mBuckets.get((int) Math.round(gcContent * 100));
    }

    public List<GCPail> getBuckets()
    {
        return mBuckets;
    }
}
