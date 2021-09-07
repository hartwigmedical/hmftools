package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvCluster;

import org.apache.commons.compress.utils.Lists;

public class CfChain
{
    public final int Id;
    public final List<CfSvData> SVs;
    public final List<CfLink> Links;

    public CfChain(int id)
    {
        Id = id;
        SVs = Lists.newArrayList();
        Links = Lists.newArrayList();
    }

    public void addSV(final CfSvData cfSvData)
    {
        SVs.add(cfSvData);
    }

    public void buildLinks()
    {
        for(final CfSvData var1 : SVs)
        {
            for(final CfSvData var2 : SVs)
            {
                if(var1 == var2)
                    continue;

                for(int se1 = SE_START; se1 <= SE_END; ++se1)
                {
                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        final int se2Index = se2;

                        if(hasAdjacentBreakends(var1, var2, se1, se2))
                        {
                            if(breakendsFace(var1, var2, se1, se2))
                            {
                                CfLink newLink = new CfLink(var1, var2, se1, se2);
                                final CfLink link1 = var1.getChainLinks()[se1];
                                final CfLink link2 = var1.getChainLinks()[se2];

                                // always keep the shortest link if there are multiple to make
                                if(link1 != null && link1.length() <= newLink.length())
                                    continue;

                                if(link2 != null && link2.length() <= newLink.length())
                                    continue;

                                Links.add(newLink);
                                var1.setChainLink(newLink, se1);
                                var2.setChainLink(newLink, se2);
                            }
                        }
                    }
                }
            }
        }
    }

    private boolean hasAdjacentBreakends(final CfSvData var1, final CfSvData var2, final int se1, final int se2)
    {
        return var1.AdjacentBreakendIds.get(se1).stream().anyMatch(x -> x == var2.BreakpointIds[se2])
                && var2.AdjacentBreakendIds.get(se2).stream().anyMatch(x -> x == var1.BreakpointIds[se1]);
    }

    private boolean breakendsFace(final CfSvData var1, final CfSvData var2, final int se1, final int se2)
    {
        if(var1.Orientations[se1] == var2.Orientations[se2])
            return false;

        if(var1.Positions[se1] < var2.Positions[se2] && var1.Orientations[se1] == NEG_ORIENT)
            return true;
        else if(var2.Positions[se1] < var1.Positions[se2] && var2.Orientations[se1] == NEG_ORIENT)
            return true;
        else
            return false;
    }

    public int getSharedSvCount(final CfSvData cfSvData)
    {
        final SvCluster svCluster = cfSvData.getSvData().getCluster();
        return (int) SVs.stream().filter(x -> x.getSvData().getCluster() == svCluster).count();
    }

    public int findMinBreakendDistance(final CfSvData cfSvData)
    {
        int minDistance = -1;

        for(final CfSvData other : SVs)
        {
            if(other == cfSvData)
                continue;

            for(int se1 = SE_START; se1 <= SE_END; ++se1)
            {
                for(int se2 = SE_START; se2 <= SE_END; ++se2)
                {
                    if(cfSvData.Chromosomes[se1].equals(other.Chromosomes[se2]))
                    {
                        int distance = abs(cfSvData.Positions[se1] - other.Positions[se2]);

                        minDistance = minDistance == -1 ? distance : min(distance, minDistance);
                    }
                }
            }
        }

        return minDistance;
    }

    public String toString()
    {
        return String.format("%d: SV(%s) links(%d)", Id, SVs.size(), Links.size()); 
    }

}
