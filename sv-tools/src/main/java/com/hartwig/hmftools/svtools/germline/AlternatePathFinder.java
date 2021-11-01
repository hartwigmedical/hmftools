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

    public static List<AlternatePath> findPaths(final List<Link> links, final List<SvData> svList)
    {
        Set<String> failed = Sets.newHashSet();
        Map<String,AlternatePath> result = Maps.newHashMap();

        TransitiveLink transitiveLink = new TransitiveLink(svList);

        for(SvData sv : svList)
        {
            if(sv.isSgl())
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = sv.breakends()[se];
                Breakend otherBreakend = sv.breakends()[switchIndex(se)];

                if(result.containsKey(otherBreakend.VcfId) && failed.contains(otherBreakend.VcfId))
                    continue;




            }
        }



        return result.values().stream().collect(Collectors.toList());
    }
    /*
            operator fun invoke(assemblyLinkStore: LinkStore, variantStore: VariantStore): Collection<AlternatePath> {

            val failed = HashSet<String>()
            val result = HashMap<String, AlternatePath>()
            val transitiveLink = TransitiveLink(assemblyLinkStore, variantStore)

            for (variant in variantStore.selectAll()) {
                if (variant.mateId != null && !result.keys.contains(variant.mateId) && !failed.contains(variant.mateId)) {
                    val links = transitiveLink.transitiveLink(variant)
                    if (links.isNotEmpty()) {
                        val alternatePath = AlternatePath(variant.vcfId, variant.mateId, links)
                        val reverseAlternatePath = AlternatePath(variant.mateId, variant.vcfId, links.map { x -> x.reverse() }.reversed())
                        result[variant.vcfId] = alternatePath
                        result[variant.mateId] = reverseAlternatePath
                        logger.debug("Found alternate mapping of $variant CIPOS:${variant.confidenceInterval} IMPRECISE:${variant.imprecise} -> ${alternatePath.pathString()}")
                    } else {
                        failed.add(variant.vcfId)
                    }
                }
            }

            return result.values
        }
     */
}
