package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.esvee.common.CommonUtils.withinLineProximity;

import java.util.List;
import java.util.Map;

public final class LineChecker
{
    public static void markLineSites(final Map<String,List<Breakend>> chrBreakendMap)
    {
        for(List<Breakend> breakends : chrBreakendMap.values())
        {
            for(int i = 0; i < breakends.size() - 1; ++i)
            {
                Breakend breakend = breakends.get(i);
                Breakend nextBreakend = breakends.get(i + 1);

                if(breakend.otherBreakend() == nextBreakend)
                    continue;

                if(!withinLineProximity(breakend.Position, nextBreakend.Position, breakend.Orient, nextBreakend.Orient))
                    continue;

                // mark if either site is line
                if(breakend.isLine() || nextBreakend.isLine())
                {
                    breakend.setLineSiteBreakend(nextBreakend);
                    nextBreakend.setLineSiteBreakend(breakend);
                }
            }
        }
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
