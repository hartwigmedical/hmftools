package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

public class AssemblyLinks
{
    public static LinkStore buildAssembledLinks(final List<SvData> svList)
    {
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

        GR_LOGGER.debug("found {} unique assemblies", assemblyBreakendMap.size());

        LinkStore assemblyLinkStore = new LinkStore();

        for(Map.Entry<String,List<Breakend>> entry : assemblyBreakendMap.entrySet())
        {
            String assembly = entry.getKey();
            List<Breakend> breakends = entry.getValue();

            if(breakends.size() < 2)
                continue;

            int linkCounter = 0;
            for(int i = 0; i < breakends.size() - 1; ++i)
            {
                Breakend breakend1 = breakends.get(i);

                for(int j = i + 1; j < breakends.size() - 1; ++j)
                {
                    Breakend breakend2 = breakends.get(j);

                    if(breakend1.sv() == breakend2.sv())
                        continue;

                    String linkId = String.format("%d-%d", assembly, linkCounter++);

                    Link newLink = Link.from(linkId, breakend1, breakend2);
                    assemblyLinkStore.addLink(breakend1, newLink);
                    assemblyLinkStore.addLink(breakend2, newLink);
                }
            }
        }

        return assemblyLinkStore;
    }
}
