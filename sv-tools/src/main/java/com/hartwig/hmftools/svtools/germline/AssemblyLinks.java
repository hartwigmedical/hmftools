package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class AssemblyLinks
{
    public static List<Link> buildAssembledLinks(final List<SvData> svList)
    {
        List<Link> links = Lists.newArrayList();

        Map<String,List<Breakend>> assemblyBreakendMap = Maps.newHashMap();

        for(SvData sv : svList)
        {
            if(sv.isSgl())
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = sv.breakends()[se];

                for(String assembly : breakend.getAssemblies())
                {
                    List<Breakend> breakends = assemblyBreakendMap.get(assembly);

                    if(breakends == null)
                    {
                        breakends = Lists.newArrayList();
                        assemblyBreakendMap.put(assembly, breakends);
                    }

                    breakends.add(breakend);
                }
            }
        }

        GM_LOGGER.debug("found {} unique assemblies", assemblyBreakendMap.size());

        for(Map.Entry<String,List<Breakend>> entry : assemblyBreakendMap.entrySet())
        {
            List<Breakend> breakends = entry.getValue();

            if(breakends.size() < 2)
                continue;

            for(int i = 0; i < breakends.size() - 1; ++i)
            {
                Breakend breakend1 = breakends.get(i);

                for(int j = i + 1; j < breakends.size() - 1; ++j)
                {
                    Breakend breakend2 = breakends.get(j);

                    if(breakend1.SvId.equals(breakend2.SvId))
                        continue;

                    links.add(Link.from(breakend1, breakend2));
                }
            }
        }

        return links;
    }
}
