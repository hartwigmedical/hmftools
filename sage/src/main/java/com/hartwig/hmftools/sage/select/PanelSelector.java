package com.hartwig.hmftools.sage.select;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.sage.select.ReadPanelStatus.MIXED;
import static com.hartwig.hmftools.sage.select.ReadPanelStatus.OUTSIDE_PANEL;
import static com.hartwig.hmftools.sage.select.ReadPanelStatus.WITHIN_PANEL;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.region.BaseRegion;

public class PanelSelector
{
    private final List<BaseRegion> mRegions;
    private int mIndex;

    public PanelSelector(final List<BaseRegion> regions)
    {
        mRegions = regions;
        mIndex = 0;
    }

    public ReadPanelStatus panelStatus(int position)
    {
        if(mRegions == null || mRegions.isEmpty())
            return OUTSIDE_PANEL;

        updateCurrent(position, position);

        BaseRegion current = current();

        return current.containsPosition(position) ? WITHIN_PANEL : OUTSIDE_PANEL;
    }

    public ReadPanelStatus panelStatus(int start, int end)
    {
        if(mRegions == null || mRegions.isEmpty())
            return OUTSIDE_PANEL;

        updateCurrent(start, end);

        BaseRegion current = current();

        if(positionsWithin(start, end, current.start(), current.end()))
            return WITHIN_PANEL;

        if(start > current.end() || end < current.start())
            return OUTSIDE_PANEL;

        return MIXED;
    }

    public List<BaseRegion> regions() { return mRegions; }

    public boolean inPanel(int start, int end)
    {
        // returns true if start and end fall within the same panel region
        if(mRegions == null || mRegions.isEmpty())
            return false;

        updateCurrent(start, end);

        BaseRegion current = current();

        if(start <= current.end() && end >= current.start())
            return true;

        return false;
    }

    private void updateCurrent(int start, int end)
    {
        if(mRegions == null || mRegions.isEmpty())
            return;

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
    }

    private BaseRegion current()
    {
        return mRegions.get(mIndex);
    }
}
