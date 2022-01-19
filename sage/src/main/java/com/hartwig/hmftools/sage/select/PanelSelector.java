package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.utils.sv.BaseRegion;

@NotThreadSafe
public class PanelSelector
{
    private final List<BaseRegion> mRegions;
    private int mIndex;

    public PanelSelector(final List<BaseRegion> regions)
    {
        mRegions = regions;
        mIndex = 0;
    }

    public boolean inPanel(int start, int end)
    {
        if(mRegions == null || mRegions.isEmpty())
            return false;

        BaseRegion current = current();
        while(mIndex > 0 && current.start() > end)
        {
            mIndex--;
            current = current();
        }

        while(mIndex < mRegions.size() - 1 && current.end() < start)
        {
            mIndex++;
            current = current();
        }

        if(start <= current.end() && end >= current.start())
            return true;

        return false;
    }

    private BaseRegion current()
    {
        return mRegions.get(mIndex);
    }

}
