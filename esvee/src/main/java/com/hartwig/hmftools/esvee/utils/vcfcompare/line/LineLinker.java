package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;

import org.jetbrains.annotations.Nullable;

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

    private static boolean nearbyBreakendsMeetLineCriteria(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
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

    private static boolean canLinkBreakendsAsLinePair(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        if(maybeInsertSite.hasLineLink() || maybeLinkedSite.hasLineLink())
            return false;

        if(maybeInsertSite == maybeLinkedSite)
            return false;

        return nearbyBreakendsMeetLineCriteria(maybeInsertSite, maybeLinkedSite);
    }

    private static @Nullable LineLink tryLinkLineBreakendPair(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite, LineLinkType linkType)
    {
        LineLink lineLink = null;

        if(canLinkBreakendsAsLinePair(maybeInsertSite, maybeLinkedSite))
        {
            lineLink = new LineLink(maybeInsertSite, maybeLinkedSite, linkType);

            maybeInsertSite.LinkedLineBreakends = lineLink;
            maybeLinkedSite.LinkedLineBreakends = lineLink;
        }

        return lineLink;
    }

    public static void linkBreakends(Map<String, List<VariantBreakend>> chrBreakendMap)
    {
        SV_LOGGER.info("Linking potential LINE breakends");

        int linkCount = 0;

        for(List<VariantBreakend> breakends : chrBreakendMap.values())
        {
            for(VariantBreakend breakend1 : breakends)
            {
                for(VariantBreakend breakend2 : breakends)
                {
                    LineLink lineLink = tryLinkLineBreakendPair(breakend1, breakend2, LineLinkType.LINKED);

                    if(lineLink != null)
                        linkCount++;
                }
            }
        }

        if(linkCount > 0)
        {
            SV_LOGGER.debug("Formed {} LINE links", linkCount);
        }
    }

    public static void inferLinks(
            Map<String, List<VariantBreakend>> sourceChrBreakendMap,
            Map<String, List<VariantBreakend>> targetChrBreakendMap
    )
    {
        SV_LOGGER.info("Inferring LINE links using poly A sites from different files");

        int linkCount = 0;

        for(String chromosome : sourceChrBreakendMap.keySet())
        {
            List<VariantBreakend> sourceBreakends = sourceChrBreakendMap.get(chromosome);
            List<VariantBreakend> targetBreakends = targetChrBreakendMap.get(chromosome);

            if(targetBreakends == null)
                continue;

            for(VariantBreakend sourceBreakend : sourceBreakends)
            {
                if(!sourceBreakend.hasPolyATail())
                    continue;

                for(VariantBreakend targetBreakend : targetBreakends)
                {
                    if(targetBreakend.hasLineLink())
                        continue;

                    LineLink lineLink = tryLinkLineBreakendPair(sourceBreakend, targetBreakend, LineLinkType.INFERRED_OTHER);

                    if(lineLink != null)
                        linkCount++;
                }
            }
        }

        if(linkCount > 0)
        {
            SV_LOGGER.debug("Inferred {} LINE links", linkCount);
        }
    }
}
