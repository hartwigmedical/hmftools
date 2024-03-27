package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegions;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

public class BlacklistLocations
{
    private final Map<String,List<BaseRegion>> mChrLocationsMap; // keyed by chromosome start
    private final boolean mIsValid;

    public BlacklistLocations(final String filename)
    {
        mChrLocationsMap = loadChrBaseRegions(filename, true);
        mIsValid = mChrLocationsMap != null;
    }

    public List<BaseRegion> getRegions(final String chromosome) { return mChrLocationsMap.get(chromosome); }
    public boolean isValid() { return mIsValid; }
    public int size() { return mChrLocationsMap.values().stream().mapToInt(x -> x.size()).sum(); }

    public boolean inBlacklistLocation(final String chromosome, final int posStart, int posEnd)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(chromosome);
        return regions != null ? regions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), posStart, posEnd)) : false;
    }

    public BaseRegion findBlacklistLocation(final String chromosome, final int position)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(chromosome);
        return regions != null ? regions.stream().filter(x -> x.containsPosition(position)).findFirst().orElse(null) : null;
    }

    public void addRegion(final String chromosome, final BaseRegion region)
    {
        List<BaseRegion> regions = mChrLocationsMap.get(chromosome);

        if(regions == null)
        {
            regions = Lists.newArrayList();
            mChrLocationsMap.put(chromosome, regions);
        }

        regions.add(region);
    }
}
