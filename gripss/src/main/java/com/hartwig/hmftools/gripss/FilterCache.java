package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.filters.FilterType.DEDUP;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.HotspotCache;

public class FilterCache
{
    private final Map<Breakend,List<FilterType>> mBreakendFilters; // empty if a breakend is not filtered
    private final Map<Breakend,List<FilterType>> mRescuedBreakendFilters; // for debug, the types of filters that a breakend was rescued from
    private final Set<Breakend> mDuplicateBreakends;

    private final List<SvData> mHotspots;
    private int mPonFiltered;

    public FilterCache()
    {
        mBreakendFilters = Maps.newHashMap();
        mRescuedBreakendFilters = Maps.newHashMap();
        mHotspots = Lists.newArrayList();
        mDuplicateBreakends = Sets.newHashSet();
        mPonFiltered = 0;
    }

    public List<SvData> getHotspots() { return mHotspots; }
    public boolean isHotspot(final SvData sv) { return mHotspots.contains(sv); }

    public Map<Breakend,List<FilterType>> getBreakendFilters() { return mBreakendFilters; }
    public Map<Breakend,List<FilterType>> getRescuedBreakendFilters() { return mRescuedBreakendFilters; }

    public List<FilterType> getBreakendFilters(final Breakend breakend) { return mBreakendFilters.get(breakend); }
    public boolean hasFilters(final Breakend breakend) { return mBreakendFilters.containsKey(breakend); }

    public boolean hasFilter(final SvData sv, final FilterType filter)
    {
        List<FilterType> filters = getBreakendFilters(sv.breakendStart());
        if(filters != null && filters.contains(filter))
            return true;

        if(sv.isSgl())
            return false;

        filters = getBreakendFilters(sv.breakendStart());
        return filters != null && filters.contains(filter);
    }

    public Set<Breakend> getDuplicateBreakends() { return mDuplicateBreakends; }

    public int ponFilteredCount() { return mPonFiltered; }

    public void checkPonFilter(final PonCache ponCache, final SvData sv)
    {
        if(ponCache.getPonCount(sv) > 0)
        {
            addBreakendFilter(sv.breakendStart(), PON);

            if(!sv.isSgl())
                addBreakendFilter(sv.breakendEnd(), PON);
            ++mPonFiltered;
        }
    }

    public void checkHotspotFilter(final HotspotCache hotspotCache, final SvData sv)
    {
        if(hotspotCache.isHotspotVariant(sv))
            mHotspots.add(sv);
    }

    public void addBreakendFilters(final Breakend breakend, final List<FilterType> filters)
    {
        mBreakendFilters.put(breakend, filters);
    }

    public void addBreakendFilter(final Breakend breakend, final FilterType filter)
    {
        List<FilterType> filters = mBreakendFilters.get(breakend);

        if(filters == null)
            mBreakendFilters.put(breakend, Lists.newArrayList(filter));
        else
            filters.add(filter);
    }

    public void updateFilters(final Set<Breakend> rescuedBreakends, final Set<Breakend> duplicateBreakends)
    {
        // add duplicate filter and remove any rescued breakends
        for(Breakend breakend : rescuedBreakends)
        {
            List<FilterType> filters = mBreakendFilters.get(breakend);

            if(filters != null)
            {
                mRescuedBreakendFilters.put(breakend, filters);
                mBreakendFilters.remove(breakend);
            }
        }

        duplicateBreakends.forEach(x -> addBreakendFilter(x, DEDUP));
        mDuplicateBreakends.addAll(duplicateBreakends);
    }

    public List<FilterType> combineSvFilters(final SvData sv)
    {
        if(sv.isSgl())
            return getBreakendFilters(sv.breakendStart());

        // take the union of the filters for an SV
        List<FilterType> combinedFilters = null;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            Breakend breakend = sv.breakends()[se];

            List<FilterType> breakendFilters = mBreakendFilters.get(breakend);

            if(breakendFilters == null)
                continue;

            if(combinedFilters == null)
            {
                combinedFilters = Lists.newArrayList(breakendFilters);
                continue;
            }

            for(FilterType filter : breakendFilters)
            {
                if(!combinedFilters.contains(filter))
                    combinedFilters.add(filter);
            }
        }

        return combinedFilters;
    }

    public void logRescuedBreakendFilters()
    {
        if(!GR_LOGGER.isDebugEnabled())
            return;

        for(Map.Entry<Breakend,List<FilterType>> entry : mRescuedBreakendFilters.entrySet())
        {
            Breakend breakend = entry.getKey();
            List<FilterType> filters = entry.getValue();
            StringJoiner sj = new StringJoiner(";");
            filters.forEach(x -> sj.add(FilterType.vcfName(x)));

            GR_LOGGER.trace("breakend({}) rescued from filters({})", breakend, sj.toString());
        }
    }

    public void clear()
    {
        mBreakendFilters.clear();
        mDuplicateBreakends.clear();
        mRescuedBreakendFilters.clear();
        mHotspots.clear();
        mPonFiltered = 0;
    }
}
