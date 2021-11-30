package com.hartwig.hmftools.gripss.links;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;

public class LinkStore
{
    private final Map<Breakend,List<Link>> mBreakendLinksMap;

    public LinkStore()
    {
        mBreakendLinksMap = Maps.newHashMap();
    }

    public Map<Breakend,List<Link>> getBreakendLinksMap() { return mBreakendLinksMap; }

    public List<Link> getBreakendLinks(final Breakend breakend) { return mBreakendLinksMap.get(breakend); }

    public static LinkStore from(final LinkStore store1, final LinkStore store2)
    {
        LinkStore newStore = new LinkStore();
        newStore.getBreakendLinksMap().putAll(store1.getBreakendLinksMap());

        for(Map.Entry<Breakend,List<Link>> entry : store2.getBreakendLinksMap().entrySet())
        {
            List<Link> links = newStore.getBreakendLinksMap().get(entry.getKey());

            if(links == null)
            {
                links = Lists.newArrayList();
                newStore.getBreakendLinksMap().put(entry.getKey(), links);
            }

            List<Link> constLinks = links;
            entry.getValue().forEach(x -> constLinks.add(x));
        }

        return newStore;
    }

    public void addLink(final Breakend breakend, final Link link)
    {
        List<Link> links = mBreakendLinksMap.get(breakend);
        if(links ==  null)
        {
            links = Lists.newArrayList();
            mBreakendLinksMap.put(breakend, links);
        }

        links.add(link);
    }

}
