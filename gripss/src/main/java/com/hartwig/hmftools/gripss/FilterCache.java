package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.filters.FilterType.DEDUP;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.HotspotCache;

public class FilterCache
{
    private final Map<Breakend, List<FilterType>> mBreakendFilters; // empty if a breakend is not filtered

    private final List<SvData> mHotspots;

    public FilterCache()
    {
        mBreakendFilters = Maps.newHashMap();
        mHotspots = Lists.newArrayList();
    }

    public void checkPonFilter(final PonCache ponCache, final SvData sv)
    {
        if(ponCache.getPonCount(sv) > 0)
        {
            addBreakendFilter(sv.breakendStart(), PON);

            if(!sv.isSgl())
                addBreakendFilter(sv.breakendEnd(), PON);
        }
    }

    public void checkHotspotFilter(final HotspotCache hotspotCache, final SvData sv)
    {
        if(hotspotCache.isHotspotVariant(sv))
            mHotspots.add(sv);
    }

    public void clear()
    {
        mBreakendFilters.clear();
        mHotspots.clear();
    }

    public List<SvData> getHotspots() { return mHotspots; }
    public boolean isHotspot(final SvData sv) { return mHotspots.contains(sv); }

    public Map<Breakend,List<FilterType>> getBreakendFilters() { return mBreakendFilters; }

    public List<FilterType> getBreakendFilters(final Breakend breakend) { return mBreakendFilters.get(breakend); }
    public boolean hasFilters(final Breakend breakend) { return mBreakendFilters.containsKey(breakend); }

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

    public List<FilterType> getSvFilters(final SvData sv)
    {
        List<FilterType> combinedFilters = Lists.newArrayList();

        for(int se = SE_START; se <= SE_END; ++se)
        {
            Breakend breakend = sv.breakends()[se];

            if(breakend == null)
                continue;

            List<FilterType> breakendFilters = mBreakendFilters.get(breakend);

            if(breakendFilters == null)
                continue;

            breakendFilters.stream().filter(x -> !combinedFilters.contains(x)).forEach(x -> combinedFilters.add(x));
        }

        return combinedFilters;
    }

    public void updateFilters(final Set<Breakend> rescuedBreakends, final Set<Breakend> duplicateBreakends)
    {
        // add duplicate filter and remove any rescued breakends
        rescuedBreakends.forEach(x -> mBreakendFilters.remove(x));
        duplicateBreakends.forEach(x -> addBreakendFilter(x, DEDUP));
    }
}
