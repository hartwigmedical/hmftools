package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.gripss.GripssConstants.MAX_DSB_DISTANCE;
import static com.hartwig.hmftools.gripss.GripssConstants.MAX_DSB_SEEK_DISTANCE;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.gripss.SvDataCache;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

public final class DsbLinkFinder
{
    // a breakend can only be in one DSB - the first & closest found

    public static LinkStore findBreaks(final SvDataCache dataCache, final LinkStore assemblyLinks, final Set<Breakend> duplicateBreakends)
    {
        LinkStore dsbLinks = new LinkStore();
        int linkId = 1;

        for(SvData sv : dataCache.getSvList())
        {
            for(Breakend breakend : sv.breakends())
            {
                if(breakend == null || dsbLinks.getBreakendLinks(breakend) != null)
                    continue;

                // ignore duplicate breakends
                if(duplicateBreakends.contains(breakend))
                    continue;

                Link[] links = findLinks(linkId++, breakend, dataCache, assemblyLinks, dsbLinks, duplicateBreakends);

                if(links == null || links.length != 2)
                    continue;

                dsbLinks.addLink(links[0].breakendStart(), links[0]);
                dsbLinks.addLink(links[1].breakendStart(), links[1]);
            }
        }

        return dsbLinks;
    }

    private static Link[] findLinks(
            int linkId, final Breakend breakend, final SvDataCache dataCache, final LinkStore assemblyLinks, final LinkStore dsbLinks,
            final Set<Breakend> duplicateBreakends)
    {
        List<Breakend> nearbyBreakends = dataCache.selectOthersNearby(breakend, MAX_DSB_DISTANCE, MAX_DSB_SEEK_DISTANCE).stream()
                .filter(x -> x.Orientation != breakend.Orientation)
                .filter(x -> !duplicateBreakends.contains(x))
                .filter(x -> dsbLinks.getBreakendLinks(x) == null)
                .collect(Collectors.toList());

        if(nearbyBreakends.size() != 1)
            return null;

        Breakend otherBreakend = nearbyBreakends.get(0);

        // check the other breakend can only make a DSB with this breakend and not others too
        List<Breakend> otherNearbyBreakends = dataCache.selectOthersNearby(otherBreakend, MAX_DSB_DISTANCE, MAX_DSB_SEEK_DISTANCE).stream()
                .filter(x -> x.Orientation != otherBreakend.Orientation)
                .filter(x -> !duplicateBreakends.contains(x))
                .collect(Collectors.toList());

        if(otherNearbyBreakends.size() != 1 || otherNearbyBreakends.get(0) != breakend)
            return null;

        List<Link> existingAssemblyLinks = assemblyLinks.getBreakendLinks(breakend);

        // ignore if already linked
        if(existingAssemblyLinks != null && existingAssemblyLinks.stream().anyMatch(x -> x.otherBreakend(breakend) == otherBreakend))
            return null;

        Link[] links = new Link[2];
        String linkStr = String.format("dsb%d", linkId);
        links[0] = Link.from(linkStr, breakend, otherBreakend);
        links[1] = Link.from(linkStr, otherBreakend, breakend);
        return links;
    }
}
