package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.purple.CopyNumberMethod;

public abstract class ExtendRegion
{
    private final CopyNumberMethod mMethod;

    ExtendRegion(final CopyNumberMethod method)
    {
        this.mMethod = method;
    }

    public List<CombinedRegion> extend(final List<CombinedRegion> regions)
    {
        for(int i = 0; i < regions.size(); i++)
        {
            CombinedRegion region = regions.get(i);
            if(region.copyNumberMethod().equals(mMethod))
            {
                extendRight(i, regions);
                i -= extendLeft(i, regions);
            }
        }

        return regions;
    }

    protected abstract void extend(final CombinedRegion target, final CombinedRegion neighbour);

    private void extendRight(int startIndex, final List<CombinedRegion> regions)
    {
        assert (startIndex < regions.size());
        final CombinedRegion target = regions.get(startIndex);
        int targetIndex = startIndex + 1;

        while(targetIndex < regions.size())
        {
            final CombinedRegion neighbour = regions.get(targetIndex);

            if(ExtendUtils.doNotExtend(target, neighbour.region()))
            {
                break;
            }

            if(!neighbour.isProcessed())
            {
                target.extend(neighbour.region());
            }
            else if(neighbour.copyNumberMethod().equals(mMethod))
            {
                extend(target, neighbour);
            }
            else
            {
                break;
            }

            regions.remove(targetIndex);
        }
    }

    private int extendLeft(int startIndex, final List<CombinedRegion> regions)
    {
        assert (startIndex < regions.size());
        final CombinedRegion target = regions.get(startIndex);

        int targetIndex = startIndex - 1;
        while(targetIndex >= 0)
        {
            final CombinedRegion neighbour = regions.get(targetIndex);
            if(ExtendUtils.doNotExtend(target, neighbour.region()))
            {
                break;
            }

            if(!neighbour.isProcessed())
            {
                target.extend(neighbour.region());
            }
            else if(neighbour.copyNumberMethod().equals(mMethod))
            {
                extend(target, neighbour);
            }
            else
            {
                break;
            }

            regions.remove(targetIndex);
            targetIndex--;
        }

        return startIndex - targetIndex - 1;
    }
}
