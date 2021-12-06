package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.Link;
import com.hartwig.hmftools.gripss.links.LinkStore;

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

    public void findDuplicateSVs(final List<AlternatePath> alternatePaths)
    {
        for(AlternatePath altPath : alternatePaths)
        {
            boolean firstIsPass = !mFilterCache.hasFilters(altPath.First);

            boolean anyInAltPathPasses = altPath.Links.stream()
                    .anyMatch(x -> !mFilterCache.hasFilters(x.breakendStart()) || !mFilterCache.hasFilters(x.breakendEnd()));

            if(altPath.Links.size() == 1)
            {
                Breakend first = altPath.First;
                Breakend second = altPath.Links.get(0).breakendEnd();

                if(!keepOriginal(first, second, firstIsPass, anyInAltPathPasses))
                {
                    GR_LOGGER.trace("breakend({}) duplicate vs other({})", first, second);
                    mDuplicateBreakends.add(first);
                }

                // no need to rescue the linked variant since is already passing (ie anyInAltPathPasses is true)
            }
            else
            {
                GR_LOGGER.trace("SV({}) duplicate vs alt-path links({})", altPath.First.sv(), altPath.pathString());

                mDuplicateBreakends.add(altPath.First);
                mDuplicateBreakends.add(altPath.Second);

                if(firstIsPass || anyInAltPathPasses)
                {
                    for(Link link : altPath.Links)
                    {
                        mRescueBreakends.add(link.breakendStart());
                        mRescueBreakends.add(link.breakendEnd());
                    }
                }
            }
        }

        // duplicates supersede rescues
        mDuplicateBreakends.forEach(x -> mRescueBreakends.remove(x));
    }

    private static boolean keepOriginal(final Breakend original, final Breakend other, boolean originalIsPass, boolean otherIsPass)
    {
        // Favour PRECISE, PASSING, then QUAL
        if(original.imprecise() != other.imprecise())
            return !original.imprecise();

        if(originalIsPass != otherIsPass)
            return originalIsPass;

        return original.Qual > other.Qual;
    }

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

            for(Breakend otherBreakend : nearbyBreakends)
            {
                if(!isDuplicateCandidate(breakend, otherBreakend))
                    continue;

                if(!keepSingle(isPass, breakend, otherBreakend, linkStore))
                {
                    GR_LOGGER.trace("breakend({}) duplicate vs other({})", breakend, otherBreakend);
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
                for(Breakend otherBreakend : nearbyBreakends)
                {
                    if(!isDuplicateCandidate(breakend, otherBreakend))
                        continue;

                    GR_LOGGER.trace("breakend({}) duplicate vs other({})", otherBreakend, breakend);
                    mSingleDuplicates.add(otherBreakend);

                    if(!otherBreakend.isSgl())
                        mSingleDuplicates.add(otherBreakend.otherBreakend());
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
