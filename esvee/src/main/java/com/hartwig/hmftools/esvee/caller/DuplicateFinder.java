package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;

public class DuplicateFinder
{
    private final FilterCache mFilterCache;
    private final SvDataCache mDataCache;
    
    private final Set<Breakend> mDuplicateBreakends;
    private final Set<Breakend> mRescueBreakends;
    private final Set<Breakend> mSingleDuplicates; // SGLs which are duplicates or duplicates of SGLs

    private static final int MAX_DEDUP_SGL_SEEK_DISTANCE = 1000;
    private static final int MAX_DEDUP_SGL_ADDITIONAL_DISTANCE = 0;

    public DuplicateFinder(final SvDataCache dataCache, final FilterCache filterCache)
    {
        mDataCache = dataCache;
        mFilterCache = filterCache;
        
        mDuplicateBreakends = Sets.newHashSet();
        mRescueBreakends = Sets.newHashSet();
        mSingleDuplicates = Sets.newHashSet();
    }

    public Set<Breakend> duplicateBreakends() { return mDuplicateBreakends; }
    public Set<Breakend> rescueBreakends() { return mRescueBreakends; }
    public Set<Breakend> duplicateSglBreakends() { return mSingleDuplicates; };

    public void findDuplicateSingles(final LinkStore linkStore)
    {
        for(SvData sv : mDataCache.getSvList())
        {
            if(!sv.isSgl())
                continue;
            
            Breakend breakend = sv.breakendStart();
            
            boolean isPass = !mFilterCache.hasFilters(breakend);

            List<Breakend> nearbyBreakends = mDataCache.selectOthersNearby(
                    breakend, MAX_DEDUP_SGL_ADDITIONAL_DISTANCE, MAX_DEDUP_SGL_SEEK_DISTANCE);

            // look through duplicate breakends in the vacinity
            // if none of them require keeping the single, then mark it as a duplicate
            boolean keepSingle = true;

            for(Breakend nearBreakend : nearbyBreakends)
            {
                if(!isDuplicateCandidate(breakend, nearBreakend))
                    continue;

                if(!keepSingle(isPass, breakend, nearBreakend, linkStore))
                {
                    SV_LOGGER.trace("breakend({}) duplicate vs other({})", breakend, nearBreakend);
                    keepSingle = false;
                    break;
                }
            }

            if(!keepSingle)
            {
                mSingleDuplicates.add(breakend);
            }
            else
            {
                for(Breakend nearBreakend : nearbyBreakends)
                {
                    if(!isDuplicateCandidate(breakend, nearBreakend))
                        continue;

                    SV_LOGGER.trace("breakend({}) duplicate vs other({})", nearBreakend, breakend);
                    mSingleDuplicates.add(nearBreakend);

                    if(!nearBreakend.isSgl())
                        mSingleDuplicates.add(nearBreakend.otherBreakend());
                }
            }
        }
    }

    private static boolean isDuplicateCandidate(final Breakend breakend, final Breakend otherBreakend)
    {
        return breakend.Orientation == otherBreakend.Orientation && (!otherBreakend.imprecise() || isExactPosition(breakend, otherBreakend));
    }

    private static boolean isExactPosition(final Breakend breakend, final Breakend otherBreakend)
    {
        return otherBreakend.Position >= breakend.minPosition() && otherBreakend.Position <= breakend.maxPosition();
    }

    private boolean keepSingle(
            boolean originalIsPass, final Breakend original, final Breakend alternative, final LinkStore linkStore)
    {
        if(linkStore.getBreakendLinks(alternative) != null)
            return false;
        
        boolean altIsPass = !mFilterCache.hasFilters(alternative);

        if(originalIsPass != altIsPass)
            return originalIsPass;

        return original.Qual > alternative.Qual;
    }

    public void clear()
    {
        mDuplicateBreakends.clear();
        mSingleDuplicates.clear();
        mRescueBreakends.clear();
    }
}
