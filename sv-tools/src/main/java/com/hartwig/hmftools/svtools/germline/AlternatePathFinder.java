package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class AlternatePathFinder
{
    public static List<AlternatePath> findPaths(final SvDataCache svDataCache, final AssemblyLinks assemblyLinks)
    {
        Set<String> failed = Sets.newHashSet();
        Map<String,AlternatePath> alternatePaths = Maps.newHashMap();

        TransitiveLinkFinder transitiveLinkFinder = new TransitiveLinkFinder(svDataCache, assemblyLinks);

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

                    // TODO - reverse links
                    AlternatePath reverseAltPath = new AlternatePath(otherBreakend.VcfId, breakend.VcfId, transLinks);
                    // val reverseAlternatePath = AlternatePath(variant.mateId, variant.vcfId, links.map { x -> x.reverse() }.reversed())

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
}
