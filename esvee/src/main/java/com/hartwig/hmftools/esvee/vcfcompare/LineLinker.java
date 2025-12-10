package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class LineLinker
{
    private static final int MIN_POLY_A_OR_T_LENGTH = 10;
    private static final String POLY_A_SEQUENCE = "A".repeat(MIN_POLY_A_OR_T_LENGTH);
    private static final String POLY_T_SEQUENCE = "T".repeat(MIN_POLY_A_OR_T_LENGTH);

    private static final int POLY_A_TO_OTHER_SITE_UPPER_DISTANCE = 40;
    private static final int POLY_A_TO_OTHER_SITE_LOWER_DISTANCE = 30;

    public static boolean hasPolyATail(final Breakend breakend)
    {
        return (breakend.Orient.isReverse() && breakend.InsertSequence.endsWith(POLY_A_SEQUENCE))
            || (breakend.Orient.isForward() && breakend.InsertSequence.startsWith(POLY_T_SEQUENCE));
    }

    private static boolean breakendPairMeetsLineCriteria(final Breakend maybeInsertSite, final Breakend maybeLinkedSite)
    {
        boolean meetsPolyACriteria = (
                maybeInsertSite.Orient.isReverse() &&
                maybeLinkedSite.Orient.isForward() &&
                maybeInsertSite.InsertSequence.endsWith(POLY_A_SEQUENCE) &&
                maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome) &&
                positionWithin(
                        maybeLinkedSite.Position,
                        maybeInsertSite.Position - POLY_A_TO_OTHER_SITE_LOWER_DISTANCE,
                        maybeInsertSite.Position + POLY_A_TO_OTHER_SITE_UPPER_DISTANCE
                )
        );

        boolean meetsPolyTCriteria = (
                maybeInsertSite.Orient.isForward() &&
                maybeLinkedSite.Orient.isReverse() &&
                maybeInsertSite.InsertSequence.startsWith(POLY_T_SEQUENCE) &&
                maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome) &&
                positionWithin(
                        maybeLinkedSite.Position,
                        maybeInsertSite.Position - POLY_A_TO_OTHER_SITE_UPPER_DISTANCE,
                        maybeInsertSite.Position + POLY_A_TO_OTHER_SITE_LOWER_DISTANCE
                )
        );

        return meetsPolyACriteria || meetsPolyTCriteria;
    }

    public static Map<String,List<Breakend>> dedupBreakends(final  Map<String,List<Breakend>> chrBreakendMap)
    {
        SV_LOGGER.debug("selecting unique variants by coords");

        Map<String, List<Breakend>> chrBreakendMapDeduped = new LinkedHashMap<>();
        int selectedCount = 0;

        for(String chromosome : chrBreakendMap.keySet())
        {
            List<Breakend> chrBreakends = chrBreakendMap.get(chromosome);

            // Group breakends by exact coords
            Map<String, List<Breakend>> coordsBreakendsMap = new LinkedHashMap<>();
            for(Breakend breakend : chrBreakends)
            {
                String coords = breakend.coordStr();

                coordsBreakendsMap.putIfAbsent(coords, new ArrayList<>());
                coordsBreakendsMap.get(coords).add(breakend);
            }

            // Dedup each breakend group
            List<Breakend> chrBreakendsDeduped = new ArrayList<>();

            for(List<Breakend> sameCoordBreakends : coordsBreakendsMap.values())
            {
                Breakend selectedBreakend = sameCoordBreakends.get(0);

                if(sameCoordBreakends.size() > 1)
                {
                    for(int nextIndex = 1; nextIndex < sameCoordBreakends.size(); nextIndex++)
                    {
                        Breakend nextBreakend = sameCoordBreakends.get(nextIndex);

                        // Prefer variant with poly A, then highest qual
                        if((nextBreakend.hasPolyATail() && !selectedBreakend.hasPolyATail())
                        || nextBreakend.sv().qual() > selectedBreakend.sv().qual())
                        {
                            selectedBreakend = nextBreakend;
                        }
                    }
                }

                chrBreakendsDeduped.add(selectedBreakend);
                selectedCount++;
            }

            // Add the deduped breakends for this chromosome to the map
            chrBreakendMapDeduped.put(chromosome, chrBreakendsDeduped);
        }

        SV_LOGGER.trace("found unique {} variants preferring poly A, then highest qual", selectedCount);

        return chrBreakendMapDeduped;
    }

    public static void linkBreakends(final Map<String,List<Breakend>> chrBreakendMap)
    {
        SV_LOGGER.info("linking breakends with LINE characteristics");

        int linkCount = 0;

        for(List<Breakend> breakends : chrBreakendMap.values())
        {
            for(Breakend maybePolyASite : breakends)
            {
                for(Breakend maybeOtherSite : breakends)
                {
                    LineLink lineLink = null;

                    if(maybePolyASite == maybeOtherSite)
                        continue;

                    if(maybePolyASite.hasLineLink() || maybeOtherSite.hasLineLink())
                        continue;

                    if(breakendPairMeetsLineCriteria(maybePolyASite, maybeOtherSite))
                    {
                        lineLink = new LineLink(maybePolyASite, maybeOtherSite, LineLinkType.LINKED);
                        maybePolyASite.setLineLink(lineLink, false);
                        maybeOtherSite.setLineLink(lineLink, false);
                    }

                    if(lineLink != null)
                        linkCount++;
                }
            }
        }

        if(linkCount > 0)
        {
            SV_LOGGER.trace("formed {} LINE links", linkCount);
        }
    }

    public static void inferLinksBetweenBreakendSets(
            final Map<String,List<Breakend>> chrMaybePolyASitesMap, final Map<String,List<Breakend>> chrMaybeOtherSitesMap,
            final LineLinkType linkType)
    {
        SV_LOGGER.debug("inferring LINE links between sets of breakends: {}", linkType);

        int linkCount = 0;

        for(String chromosome : chrMaybePolyASitesMap.keySet())
        {
            List<Breakend> maybePolyASitesList = chrMaybePolyASitesMap.get(chromosome);
            List<Breakend> maybeOtherSitesList = chrMaybeOtherSitesMap.get(chromosome);

            if(maybeOtherSitesList == null)
                continue;

            for(Breakend maybePolyASite : maybePolyASitesList)
            {
                if(!maybePolyASite.hasPolyATail())
                    continue;

                for(Breakend maybeOtherSite : maybeOtherSitesList)
                {
                    LineLink lineLink = null;

                    if(maybeOtherSite.hasLineLink())
                        continue;

                    if(breakendPairMeetsLineCriteria(maybePolyASite, maybeOtherSite))
                    {
                        lineLink = new LineLink(maybePolyASite, maybeOtherSite, linkType);
                        maybePolyASite.setLineLink(lineLink, true);
                    }

                    if(lineLink != null)
                        linkCount++;
                }
            }
        }

        if(linkCount > 0)
        {
            SV_LOGGER.trace("inferred {} LINE links", linkCount);
        }
    }
}
