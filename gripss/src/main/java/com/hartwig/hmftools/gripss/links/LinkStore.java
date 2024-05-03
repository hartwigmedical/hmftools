package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.LINKED_BY_DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

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

    public void addLinks(final String linkId, final Breakend breakend1, final Breakend breakend2)
    {
        addLinks(linkId, breakend1, breakend2, true);
    }

    public void addLinks(final String linkId, final Breakend breakend1, final Breakend breakend2, boolean allowDuplicates)
    {
        if(!allowDuplicates)
        {
            List<Link> existingLinks = getBreakendLinks(breakend1);

            if(existingLinks != null && existingLinks.stream().anyMatch(x -> x.otherBreakend(breakend1) == breakend2))
                return;
        }

        addLink(breakend1, Link.from(linkId, breakend1, breakend2));
        addLink(breakend2, Link.from(linkId, breakend2, breakend1));
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

    public String getBreakendLinksStr(final Breakend breakend)
    {
        List<Link> links = mBreakendLinksMap.get(breakend);

        if(links == null || links.isEmpty())
            return "";

        StringJoiner sj = new StringJoiner(LINKED_BY_DELIM);
        links.forEach(x -> sj.add(x.Id));
        return sj.toString();
    }

    public void clear()
    {
        mBreakendLinksMap.clear();
    }
}
