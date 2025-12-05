package com.hartwig.hmftools.cobalt.calculations;

import java.util.ArrayList;
import java.util.List;

class GCPailsList
{
    private final List<GCPail> mBuckets = new ArrayList<>(101);

    GCPailsList()
    {
        for(int i = 0; i < 101; i++)
        {
            mBuckets.add(new GCPail(i));
        }
    }

    GCPail getGCPail(double gcContent)
    {
        return mBuckets.get(GCPail.bucketIndex(gcContent));
    }

    List<GCPail> getBuckets()
    {
        return mBuckets;
    }
}
