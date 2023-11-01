package com.hartwig.hmftools.common.utils.sv;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public final class SvCommonUtils
{
    public static final byte POS_ORIENT = 1;
    public static final byte NEG_ORIENT = -1;

    public static void mergeChrBaseRegionOverlaps(final List<ChrBaseRegion> regions, boolean checkSorted)
    {
        if(checkSorted)
            Collections.sort(regions);

        // first a quick check for overlaps
        boolean hasOverlaps = false;

        for(int i = 0; i < regions.size() - 1; ++i)
        {
            ChrBaseRegion region = regions.get(i);
            ChrBaseRegion next = regions.get(i + 1);

            if(region.overlaps(next))
            {
                hasOverlaps = true;
                break;
            }
        }

        if(!hasOverlaps)
            return;

        int i = 0;
        while(i < regions.size() - 1)
        {
            ChrBaseRegion region = regions.get(i);

            int j = i + 1;
            while(j < regions.size())
            {
                ChrBaseRegion next = regions.get(i + 1);

                if(region.overlaps(next))
                {
                    region.setEnd(next.end());
                    regions.remove(j);
                }
                else
                {
                    ++j;
                }
            }

            ++i;
        }
    }
}
