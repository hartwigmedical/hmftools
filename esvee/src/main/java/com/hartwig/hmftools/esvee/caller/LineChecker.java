package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.esvee.common.CommonUtils.withinLineProximity;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.assembly.types.Junction;

public final class LineChecker
{
    public static void markLineSites(final Map<String,List<Breakend>> chrBreakendMap)
    {
        int maxProximateLineDistance = max(LINE_INDEL_MAX_OVERLAP, LINE_INDEL_MAX_GAP);

        for(List<Breakend> breakends : chrBreakendMap.values())
        {
            for(int i = 0; i < breakends.size() - 1; ++i)
            {
                Breakend breakend = breakends.get(i);

                // first search for likely remote line source sites
                for(int j = i + 1; j < breakends.size(); ++j)
                {
                    Breakend nextBreakend = breakends.get(j);

                    if(abs(breakend.Position - nextBreakend.Position) > maxProximateLineDistance)
                        break;

                    if(breakend.otherBreakend() == nextBreakend)
                        continue;

                    if(!breakend.isLine() && !nextBreakend.isLine())
                        continue;

                    // mark remote source sites
                    if(!breakend.sv().isLineSite() && !closeToOriginalJunction(breakend, maxProximateLineDistance))
                    {
                        breakend.markLineRemoteSourceSite();
                    }
                    else if(!nextBreakend.sv().isLineSite() && !closeToOriginalJunction(nextBreakend, maxProximateLineDistance))
                    {
                        nextBreakend.markLineRemoteSourceSite();
                    }
                }

                // otherwise check for pairs of line and the other insertion breakend
                Breakend nextBreakend = breakends.get(i + 1);

                if(breakend.otherBreakend() == nextBreakend)
                    continue;

                if(!breakend.isLine() && !nextBreakend.isLine())
                    continue;

                if(breakend.isLineRemoteSourceSite() || nextBreakend.isLineRemoteSourceSite())
                    continue;

                if(!withinLineProximity(breakend.Position, nextBreakend.Position, breakend.Orient, nextBreakend.Orient))
                    continue;

                // mark if either site is line
                breakend.setLineSiteBreakend(nextBreakend);
                nextBreakend.setLineSiteBreakend(breakend);
            }
        }
    }

    private static boolean closeToOriginalJunction(final Breakend breakend, int maxProximateLineDistance)
    {
        List<Junction> originalJunctions = breakend.sv().originalJunctions();

        for(Junction junction : originalJunctions)
        {
            if(breakend.Chromosome.equals(junction.Chromosome) && breakend.Orient == junction.Orient
            && abs(breakend.Position - junction.Position) <= maxProximateLineDistance)
            {
                return true;
            }
        }

        return false;
    }

    public static void adjustLineSites(final Variant var)
    {
        if(!var.isLineSite())
            return;

        boolean isGermline = var.isGermline();
        boolean isFiltered = var.isFiltered();

        boolean hasLineGermline = false;
        boolean hasLinePass = false;

        for(Breakend breakend : var.breakends())
        {
            if(breakend != null && breakend.lineSiteBreakend() != null)
            {
                hasLineGermline |= breakend.lineSiteBreakend().sv().isGermline();
                hasLinePass |= breakend.lineSiteBreakend().sv().isPass();
            }
        }

        if(!isGermline && hasLineGermline)
            var.markGermline();

        if(isFiltered && hasLinePass)
            var.filters().clear();
    }
}
