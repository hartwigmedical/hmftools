package com.hartwig.hmftools.esvee.vcfcompare.line;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.vcfcompare.match.VariantBreakend;

public class LineLinker
{
    private static final int MIN_POLY_A_OR_T_LENGTH = 10;
    private static final String POLY_A_SEQUENCE = "A".repeat(MIN_POLY_A_OR_T_LENGTH);
    private static final String POLY_T_SEQUENCE = "T".repeat(MIN_POLY_A_OR_T_LENGTH);

    private static final int POLY_A_TO_OTHER_SITE_UPPER_DISTANCE = 40;
    private static final int POLY_A_TO_OTHER_SITE_LOWER_DISTANCE = 30;

    public static boolean hasPolyATail(VariantBreakend breakend)
    {
        return
                (breakend.Orientation == -1 && breakend.InsertSequence.endsWith(POLY_A_SEQUENCE)) ||
                (breakend.Orientation == 1 && breakend.InsertSequence.startsWith(POLY_T_SEQUENCE));
    }

    private static boolean breakendPairMeetsLineCriteria(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        boolean meetsPolyACriteria = (
                maybeInsertSite.Orientation == -1 &&
                maybeLinkedSite.Orientation == 1 &&
                maybeInsertSite.InsertSequence.endsWith(POLY_A_SEQUENCE) &&
                maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome) &&
                positionWithin(
                        maybeLinkedSite.Position,
                        maybeInsertSite.Position - POLY_A_TO_OTHER_SITE_LOWER_DISTANCE,
                        maybeInsertSite.Position + POLY_A_TO_OTHER_SITE_UPPER_DISTANCE
                )
        );

        boolean meetsPolyTCriteria = (
                maybeInsertSite.Orientation == 1 &&
                maybeLinkedSite.Orientation == -1 &&
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

    public static Map<String, List<VariantBreakend>> dedupBreakends(Map<String, List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.debug("selecting unique variants by coords");

        Map<String, List<VariantBreakend>> chrBreakendMapDeduped = new LinkedHashMap<>();
        int selectedCount = 0;

        for(String chromosome : chrBreakendMap.keySet())
        {
            List<VariantBreakend> chrBreakends = chrBreakendMap.get(chromosome);

            // Group breakends by exact coords
            Map<String, List<VariantBreakend>> coordsBreakendsMap = new LinkedHashMap<>();
            for(VariantBreakend breakend : chrBreakends)
            {
                String coords = breakend.coordStr();

                coordsBreakendsMap.putIfAbsent(coords, new ArrayList<>());
                coordsBreakendsMap.get(coords).add(breakend);
            }

            // Dedup each breakend group
            List<VariantBreakend> chrBreakendsDeduped = new ArrayList<>();

            for(List<VariantBreakend> sameCoordBreakends : coordsBreakendsMap.values())
            {
                VariantBreakend selectedBreakend = sameCoordBreakends.get(0);

                if(sameCoordBreakends.size() > 1)
                {
                    for(int nextIndex = 1; nextIndex < sameCoordBreakends.size(); nextIndex++)
                    {
                        VariantBreakend nextBreakend = sameCoordBreakends.get(nextIndex);

                        // Prefer variant with poly A, then highest qual
                        if((nextBreakend.hasPolyATail() && !selectedBreakend.hasPolyATail()) ||
                                nextBreakend.qual() > selectedBreakend.qual())
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

    public static void linkBreakends(Map<String, List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.info("linking breakends with LINE characteristics");

        int linkCount = 0;

        for(List<VariantBreakend> breakends : chrBreakendMap.values())
        {
            for(VariantBreakend maybePolyASite : breakends)
            {
                for(VariantBreakend maybeOtherSite : breakends)
                {
                    LineLink lineLink = null;

                    if(maybePolyASite == maybeOtherSite)
                        continue;

                    if(maybePolyASite.hasLineLink() || maybeOtherSite.hasLineLink())
                        continue;

                    if(breakendPairMeetsLineCriteria(maybePolyASite, maybeOtherSite))
                    {
                        lineLink = new LineLink(maybePolyASite, maybeOtherSite, LineLinkType.LINKED);
                        maybePolyASite.LinkedLineBreakends = lineLink;
                        maybeOtherSite.LinkedLineBreakends = lineLink;
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
            final Map<String,List<VariantBreakend>> chrMaybePolyASitesMap, final Map<String,List<VariantBreakend>> chrMaybeOtherSitesMap,
            final LineLinkType linkType)
    {
        SV_LOGGER.debug("inferring LINE links between sets of breakends: {}", linkType);

        int linkCount = 0;

        for(String chromosome : chrMaybePolyASitesMap.keySet())
        {
            List<VariantBreakend> maybePolyASitesList = chrMaybePolyASitesMap.get(chromosome);
            List<VariantBreakend> maybeOtherSitesList = chrMaybeOtherSitesMap.get(chromosome);

            if(maybeOtherSitesList == null)
                continue;

            for(VariantBreakend maybePolyASite : maybePolyASitesList)
            {
                if(!maybePolyASite.hasPolyATail())
                    continue;

                for(VariantBreakend maybeOtherSite : maybeOtherSitesList)
                {
                    LineLink lineLink = null;

                    if(maybeOtherSite.hasLineLink())
                        continue;

                    if(breakendPairMeetsLineCriteria(maybePolyASite, maybeOtherSite))
                    {
                        lineLink = new LineLink(maybePolyASite, maybeOtherSite, linkType);
                        maybePolyASite.InferredLinkedLineBreakends = lineLink;
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
