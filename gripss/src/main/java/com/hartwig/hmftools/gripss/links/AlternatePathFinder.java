package com.hartwig.hmftools.gripss.links;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.SvDataCache;

public class AlternatePathFinder
{
    public static List<AlternatePath> findPaths(final SvDataCache svDataCache, final LinkStore assemblyLinkStore)
    {
        Map<String,AlternatePath> alternatePaths = Maps.newHashMap();

        TransitiveLinkFinder transitiveLinkFinder = new TransitiveLinkFinder(svDataCache, assemblyLinkStore);

        for(SvData sv : svDataCache.getSvList())
        {
            if(sv.isSgl())
                continue;

            Breakend breakend = sv.breakendStart();
            Breakend otherBreakend = sv.breakendEnd();

            List<Link> transLinks = transitiveLinkFinder.findTransitiveLinks(breakend);

            if(!transLinks.isEmpty())
            {
                AlternatePath altPath = new AlternatePath(breakend, otherBreakend, transLinks);

                List<Link> reversedLinks = Lists.newArrayList();
                transLinks.forEach(x -> reversedLinks.add(0, x.reverse()));
                AlternatePath reverseAltPath = new AlternatePath(otherBreakend, breakend, reversedLinks);

                alternatePaths.put(breakend.VcfId, altPath);
                alternatePaths.put(otherBreakend.VcfId, reverseAltPath);
            }
        }

        return alternatePaths.values().stream().collect(Collectors.toList());
    }

    public static Map<Breakend,String> createPathMap(final List<AlternatePath> alternatePaths)
    {
        Map<Breakend,String> idPathMap = Maps.newHashMap();

        for(AlternatePath altPath : alternatePaths)
        {
            idPathMap.put(altPath.First, altPath.pathString());
        }

        return idPathMap;
    }

    public static LinkStore createLinkStore(final List<AlternatePath> alternatePaths)
    {
        LinkStore linkStore = new LinkStore();

        for(AlternatePath altPath : alternatePaths)
        {
            List<Link> transLinks = altPath.transitiveLinks();
            transLinks.forEach(x -> linkStore.addLink(altPath.First, x));
        }

        return linkStore;
    }

}
