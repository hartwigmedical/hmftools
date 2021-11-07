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
    // TODO: store by breakend instead?
    private final Map<String,List<Link>> mBreakendLinksMap;

    public AssemblyLinks()
    {
        mBreakendLinksMap = Maps.newHashMap();
    }

    public Map<String,List<Link>> getBreakendLinksMap() { return mBreakendLinksMap; }

    public List<Link> getBreakendLinks(final Breakend breakend) { return mBreakendLinksMap.get(breakend.VcfId); }

    public void buildAssembledLinks(final List<SvData> svList)
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

        GM_LOGGER.debug("found {} unique assemblies", assemblyBreakendMap.size());

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
                    addLink(breakend1, newLink);
                    addLink(breakend2, newLink);
                }
            }
        }
    }

    private void addLink(final Breakend breakend, final Link link)
    {
        List<Link> links = mBreakendLinksMap.get(breakend.VcfId);
        if(links ==  null)
        {
            links = Lists.newArrayList();
            mBreakendLinksMap.put(breakend.VcfId, links);
        }

        links.add(link);
    }
}
