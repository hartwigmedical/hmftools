package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;

public class LineLinker
{
    private static final int MIN_POLY_A_OR_T_LENGTH = 10;
    private static final String POLY_A_SEQUENCE = "A".repeat(MIN_POLY_A_OR_T_LENGTH);
    private static final String POLY_T_SEQUENCE = "T".repeat(MIN_POLY_A_OR_T_LENGTH);

    private static final int OTHER_SITE_UPPER_BOUND = 30;
    private static final int OTHER_SITE_LOWER_BOUND = 10;

    public static boolean isLineInsertionSite(VariantBreakend breakend)
    {
        return
                (breakend.Orientation == -1 && breakend.InsertSequence.endsWith(POLY_A_SEQUENCE)) ||
                (breakend.Orientation == 1 && breakend.InsertSequence.startsWith(POLY_T_SEQUENCE));
    }

    private static boolean nearbyBreakendsMeetLineCriteria(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        return
                // Insert site with forward orientation with trailing poly A
                (maybeInsertSite.Orientation == -1 &&
                        maybeLinkedSite.Orientation == 1 &&
                        maybeInsertSite.InsertSequence.endsWith(POLY_A_SEQUENCE) &&
                        positionWithin(
                                maybeLinkedSite.Position,
                                maybeInsertSite.Position - OTHER_SITE_LOWER_BOUND,
                                maybeInsertSite.Position + OTHER_SITE_UPPER_BOUND
                        )
                ) ||

                // Insert site with reverse complement orientation with leading poly T
                (maybeInsertSite.Orientation == 1 &&
                        maybeLinkedSite.Orientation == -1 &&
                        maybeInsertSite.InsertSequence.startsWith(POLY_T_SEQUENCE) &&
                        positionWithin(
                                maybeLinkedSite.Position,
                                maybeInsertSite.Position - OTHER_SITE_UPPER_BOUND,
                                maybeInsertSite.Position + OTHER_SITE_LOWER_BOUND
                        )
                );
    }

    private static boolean tryLinkSgls(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        if(maybeInsertSite.hasLineLink() || maybeLinkedSite.hasLineLink())
            return false;

        if(!(maybeInsertSite.isSingle() && maybeLinkedSite.isSingle() && maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome)))
            return false;

        if(nearbyBreakendsMeetLineCriteria(maybeInsertSite, maybeLinkedSite))
        {
            maybeInsertSite.LinkedLineBreakend = maybeLinkedSite;
            maybeLinkedSite.LinkedLineBreakend = maybeInsertSite;
            return true;
        }

        return false;
    }

    private static boolean tryLinkTranslocations(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        if(maybeInsertSite.hasLineLink() || maybeLinkedSite.hasLineLink())
            return false;

        if(!(maybeLinkedSite.isTranslocation() && maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome)))
            return false;

        if(nearbyBreakendsMeetLineCriteria(maybeInsertSite, maybeLinkedSite))
        {
            maybeInsertSite.LinkedLineBreakend = maybeLinkedSite;
            maybeLinkedSite.LinkedLineBreakend = maybeInsertSite;
            return true;
        }

        return false;
    }

    public static void tryLinkLineBreakends(Map<String, List<VariantBreakend>> chrBreakendMap)
    {
        for(List<VariantBreakend> breakends : chrBreakendMap.values())
        {
            for(VariantBreakend breakend1 : breakends)
            {
                for(VariantBreakend breakend2 : breakends)
                {
                    tryLinkSgls(breakend1, breakend2);
                    tryLinkTranslocations(breakend1, breakend2);
                }
            }
        }
    }
}
