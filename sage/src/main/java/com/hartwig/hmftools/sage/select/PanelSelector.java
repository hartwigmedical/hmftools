package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

@NotThreadSafe
public class PanelSelector
{
    private final List<ChrBaseRegion> mRegions;
    private int mIndex;

    public PanelSelector(final List<ChrBaseRegion> regions)
    {
        mRegions = regions;
        mIndex = 0;
    }

    public boolean inPanel(int start, int end)
    {
        if(mRegions.isEmpty())
            return false;

        ChrBaseRegion current = current();
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

    private ChrBaseRegion current()
    {
        return mRegions.get(mIndex);
    }

}
