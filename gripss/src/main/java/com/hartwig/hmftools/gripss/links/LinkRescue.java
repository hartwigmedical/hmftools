package com.hartwig.hmftools.gripss.links;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.gripss.GripssConstants.SHORT_RESCUE_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.DEDUP;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.gripss.FilterCache;
import com.hartwig.hmftools.gripss.common.Breakend;

public class LinkRescue
{
    private final Map<Breakend,String> mRescueInfo;

    private static final String TYPE_LINK = "link";
    private static final String TYPE_LINE_DSB = "line_dsb";

    public LinkRescue()
    {
        mRescueInfo = Maps.newHashMap();
    }

    public Map<Breakend,String> getRescueInfo() { return mRescueInfo; }

    public void findRescuedBreakends(final LinkStore linkStore, final FilterCache filterCache, boolean rescueShortSVs)
    {
        for(Breakend breakend : linkStore.getBreakendLinksMap().keySet())
        {
            if(mRescueInfo.containsKey(breakend))
                continue;

            if(!filterCache.hasFilters(breakend))
                continue;

            if(!isRescueCandidate(breakend, filterCache, rescueShortSVs))
                continue;

            // this is a candidate for link-based rescue

            Set<Breakend> linkedBreakends = findLinkVariants(breakend, linkStore);

            if(linkedBreakends.isEmpty())
                continue;

            Breakend passingBreakend = linkedBreakends.stream()
                    .filter(x -> !filterCache.hasFilters(x))
                    .filter(x -> isRescueCandidate(x, filterCache, rescueShortSVs))
                    .findFirst().orElse(null);

            if(passingBreakend != null)
            {
                for(Breakend linkedBreakend : linkedBreakends)
                {
                    if(!isRescueCandidate(linkedBreakend, filterCache, rescueShortSVs))
                        continue;

                    addBreakend(linkedBreakend, passingBreakend, TYPE_LINK);
                }
            }
        }
    }

    private void addBreakend(final Breakend breakend, final Breakend rescuingBreakend, final String type)
    {
        String rescueInfo = String.format("%s_%s", rescuingBreakend.VcfId, type);
        mRescueInfo.put(breakend, rescueInfo);

        Breakend otherBreakend = breakend.otherBreakend();

        if(otherBreakend != null)
            mRescueInfo.put(otherBreakend, rescueInfo);
    }

    private static boolean isRescueCandidate(final Breakend breakend, final FilterCache filterCache, boolean rescueShortSVs)
    {
        if(filterCache.hasFilter(breakend.sv(), DEDUP))
            return false;

        if(!rescueShortSVs && tooShortToRescue(breakend.type(), breakend.sv().length()))
            return false;

        return true;
    }

    private static Set<Breakend> findLinkVariants(final Breakend breakend, final LinkStore linkStore)
    {
        Set<Breakend> linkedBreakends = Sets.newHashSet();
        findLinkVariants(breakend, linkStore, linkedBreakends);
        return linkedBreakends;
    }

    private static void findLinkVariants(final Breakend breakend, final LinkStore linkStore, final Set<Breakend> linkedBreakends)
    {
        if(linkedBreakends.contains(breakend))
            return;

        linkedBreakends.add(breakend);

        // call recursively on all other breakends that this one is linked to, and their SV's other breakend
        List<Link> links = linkStore.getBreakendLinks(breakend);

        if(links == null)
            return;

        for(Link link : links)
        {
            Breakend linkedBreakend = link.otherBreakend(breakend);
            findLinkVariants(linkedBreakend, linkStore, linkedBreakends);

            if(linkedBreakend.otherBreakend() != null)
                findLinkVariants(linkedBreakend.otherBreakend(), linkStore, linkedBreakends);
        }
    }

    public void findRescuedDsbLineInsertions(final LinkStore linkStore, final FilterCache filterCache, double minQual)
    {
        for(Breakend breakend : linkStore.getBreakendLinksMap().keySet())
        {
            if(mRescueInfo.containsKey(breakend))
                continue;

            if(!filterCache.hasFilters(breakend))
                continue;

            if(!isLineRescueCandidate(breakend, filterCache))
                continue;

            List<Link> links = linkStore.getBreakendLinks(breakend);

            for(Link link : links)
            {
                Breakend otherBreakend = link.otherBreakend(breakend);

                if(!breakend.IsLineInsertion && !otherBreakend.IsLineInsertion)
                    continue;

                double combinedQual = breakend.Qual + otherBreakend.Qual;

                if(combinedQual < minQual)
                    continue;

                if(!isLineRescueCandidate(otherBreakend, filterCache))
                    continue;

                addBreakend(breakend, otherBreakend, TYPE_LINE_DSB);
                addBreakend(otherBreakend, breakend, TYPE_LINE_DSB);
            }
        }
    }

    public static boolean tooShortToRescue(final StructuralVariantType type, int svLength)
    {
        return (type == DEL || type == DUP || type == INS) && svLength < SHORT_RESCUE_LENGTH;
    }

    private static boolean isLineRescueCandidate(final Breakend breakend, final FilterCache filterCache)
    {
        if(filterCache.hasFilter(breakend.sv(), PON))
            return false;

        if(filterCache.hasFilter(breakend.sv(), DEDUP))
            return false;

        if(tooShortToRescue(breakend.type(), breakend.sv().length()))
            return false;

        return true;
    }
}

