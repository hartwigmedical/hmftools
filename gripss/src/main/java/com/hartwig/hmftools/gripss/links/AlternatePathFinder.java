package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.Breakend;
import com.hartwig.hmftools.gripss.SvData;
import com.hartwig.hmftools.gripss.SvDataCache;

import org.apache.commons.compress.utils.Lists;

public class AlternatePathFinder
{
    public static List<AlternatePath> findPaths(final SvDataCache svDataCache, final LinkStore assemblyLinkStore)
    {
        Set<String> failed = Sets.newHashSet();
        Map<String,AlternatePath> alternatePaths = Maps.newHashMap();

        TransitiveLinkFinder transitiveLinkFinder = new TransitiveLinkFinder(svDataCache, assemblyLinkStore);

        for(SvData sv : svDataCache.getSvList())
        {
            if(sv.isSgl())
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = sv.breakends()[se];
                Breakend otherBreakend = sv.breakends()[switchIndex(se)];

                if(alternatePaths.containsKey(otherBreakend.VcfId) && failed.contains(otherBreakend.VcfId))
                    continue;

                List<Link> transLinks = transitiveLinkFinder.findTransitiveLinks(breakend);

                if(!transLinks.isEmpty())
                {
                    AlternatePath altPath = new AlternatePath(breakend.VcfId, otherBreakend.VcfId, transLinks);

                    List<Link> reversedLinks = Lists.newArrayList();
                    transLinks.forEach(x -> reversedLinks.add(0, x.reverse()));
                    AlternatePath reverseAltPath = new AlternatePath(otherBreakend.VcfId, breakend.VcfId, reversedLinks);

                    alternatePaths.put(breakend.VcfId, altPath);
                    alternatePaths.put(otherBreakend.VcfId, reverseAltPath);

                    // result[variant.vcfId] = alternatePath
                    // result[variant.mateId] = reverseAlternatePath
                    // logger.debug("Found alternate mapping of $variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${alternatePath.pathString()}")
                }
                else
                {
                    failed.add(breakend.VcfId);
                }
            }
        }

        return alternatePaths.values().stream().collect(Collectors.toList());
    }

    public static Map<String,String> createIdToPathMap(final List<AlternatePath> alternatePaths)
    {
        Map<String,String> idPathMap = Maps.newHashMap();

        for(AlternatePath altPath : alternatePaths)
        {
            idPathMap.put(altPath.VcfId, altPath.pathString());
        }

        return idPathMap;
    }
}
