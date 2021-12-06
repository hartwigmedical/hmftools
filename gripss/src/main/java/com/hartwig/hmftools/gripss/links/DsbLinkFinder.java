package com.hartwig.hmftools.gripss.links;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.SvDataCache;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

public final class DsbLinkFinder
{
    private static final int MAX_DSB_SEEK_DISTANCE = 1000;
    private static final int MAX_DSB_DISTANCE = 30;

    public static LinkStore findBreaks(final SvDataCache dataCache, final LinkStore assemblyLinks, final Set<Breakend> duplicateBreakends)
    {
        LinkStore dsbLinks = new LinkStore();
        int linkId = 1;

        Set<Breakend> linkedBreakends = Sets.newHashSet();

        for(SvData sv : dataCache.getSvList())
        {
            for(Breakend breakend : sv.breakends())
            {
                if(breakend == null || linkedBreakends.contains(breakend))
                    continue;

                List<Link> links = findLinks(linkId++, breakend, dataCache, assemblyLinks, duplicateBreakends);

                if(links == null)
                    continue;

                dsbLinks.addLink(links.get(0).breakendStart(), links.get(0));
                dsbLinks.addLink(links.get(1).breakendStart(), links.get(1));
                linkedBreakends.add(links.get(1).breakendStart()); // to prevent evaluating the other breakend again
            }
        }

        return dsbLinks;
    }

    private static List<Link> findLinks(
            int linkId, final Breakend breakend, final SvDataCache dataCache, final LinkStore assemblyLinks,
            final Set<Breakend> duplicateBreakends)
    {
        List<Breakend> nearbyBreakends = dataCache.selectOthersNearby(breakend, MAX_DSB_DISTANCE, MAX_DSB_SEEK_DISTANCE).stream()
                .filter(x -> x.Orientation != breakend.Orientation)
                .filter(x -> !duplicateBreakends.contains(x))
                .collect(Collectors.toList());

        if(nearbyBreakends.size() != 1)
            return null;

        Breakend otherBreakend = nearbyBreakends.get(0);

        // check the other breakend can only make a DSB with this breakend and not others too
        List<Breakend> otherNearbyBreakends = dataCache.selectOthersNearby(otherBreakend, MAX_DSB_DISTANCE, MAX_DSB_SEEK_DISTANCE).stream()
                .filter(x -> x.Orientation != otherBreakend.Orientation)
                .filter(x -> !duplicateBreakends.contains(x))
                .collect(Collectors.toList());

        if(otherNearbyBreakends.size() != 1)
            return null;

        List<Link> existingLinks = assemblyLinks.getBreakendLinks(breakend);

        // ignore if already linked
        if(existingLinks != null && existingLinks.stream().anyMatch(x -> x.otherBreakend(breakend) == otherBreakend))
            return null;

        List<Link> links = Lists.newArrayList();
        String linkStr = String.format("dsb%d", linkId);
        links.add(Link.from(linkStr, breakend, otherBreakend));
        links.add(Link.from(linkStr, otherBreakend, breakend));

        return links;
    }
}
