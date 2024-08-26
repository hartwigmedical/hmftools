package com.hartwig.hmftools.esvee.utils.vcfcompare.line;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.utils.vcfcompare.common.VariantBreakend;

public class LineLinker
{
    private final Map<String, List<VariantBreakend>> mChrBreakendMap;
    private final List<LineLink> mLinkedBreakends;
    private int mLinkCount = 0;

    private static final int MIN_POLY_A_OR_T_LENGTH = 10;
    private static final String POLY_A_SEQUENCE = "A".repeat(MIN_POLY_A_OR_T_LENGTH);
    private static final String POLY_T_SEQUENCE = "T".repeat(MIN_POLY_A_OR_T_LENGTH);

    private static final int OTHER_SITE_UPPER_BOUND = 30;
    private static final int OTHER_SITE_LOWER_BOUND = 10;

    public LineLinker(Map<String, List<VariantBreakend>> chrBreakendMap)
    {
        mChrBreakendMap = chrBreakendMap;
        mLinkedBreakends = new ArrayList<>();
    }

    public static boolean hasPolyATail(VariantBreakend breakend)
    {
        return
                (breakend.Orientation == -1 && breakend.InsertSequence.endsWith(POLY_A_SEQUENCE)) ||
                (breakend.Orientation == 1 && breakend.InsertSequence.startsWith(POLY_T_SEQUENCE));
    }

    private static boolean nearbyBreakendsMeetLineCriteria(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        boolean insertSiteHasPolyA = (maybeInsertSite.Orientation == -1 &&
                maybeLinkedSite.Orientation == 1 &&
                maybeInsertSite.InsertSequence.endsWith(POLY_A_SEQUENCE) &&
                maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome) &&
                positionWithin(
                        maybeLinkedSite.Position,
                        maybeInsertSite.Position - OTHER_SITE_LOWER_BOUND,
                        maybeInsertSite.Position + OTHER_SITE_UPPER_BOUND
                )
        );

        boolean insertSiteHasPolyT = (maybeInsertSite.Orientation == 1 &&
                maybeLinkedSite.Orientation == -1 &&
                maybeInsertSite.InsertSequence.startsWith(POLY_T_SEQUENCE) &&
                maybeInsertSite.Chromosome.equals(maybeLinkedSite.Chromosome) &&
                positionWithin(
                        maybeLinkedSite.Position,
                        maybeInsertSite.Position - OTHER_SITE_UPPER_BOUND,
                        maybeInsertSite.Position + OTHER_SITE_LOWER_BOUND
                )
        );

        return insertSiteHasPolyA || insertSiteHasPolyT;
    }

    private LineLink tryLinkLineBreakendPair(VariantBreakend maybeInsertSite, VariantBreakend maybeLinkedSite)
    {
        if(maybeInsertSite.hasLineLink() || maybeLinkedSite.hasLineLink())
            return null;

        if(maybeInsertSite == maybeLinkedSite)
            return null;

        if(!nearbyBreakendsMeetLineCriteria(maybeInsertSite, maybeLinkedSite))
            return null;

        LineLink lineLink = new LineLink(maybeInsertSite, maybeLinkedSite);

        maybeInsertSite.LinkedLineBreakends = lineLink;
        maybeLinkedSite.LinkedLineBreakends = lineLink;
        mLinkCount++;

        return lineLink;
    }

    public void tryLinkLineBreakends()
    {
        SV_LOGGER.info("Linking potential LINE breakends");

        for(List<VariantBreakend> breakends : mChrBreakendMap.values())
        {
            for(VariantBreakend breakend1 : breakends)
            {
                for(VariantBreakend breakend2 : breakends)
                {
                    LineLink lineLink = tryLinkLineBreakendPair(breakend1, breakend2);

                    if(lineLink != null)
                    {
                        mLinkedBreakends.add(lineLink);
                    }
                }
            }
        }

        if(mLinkCount > 0)
        {
            SV_LOGGER.debug("Formed {} LINE links", mLinkCount);
        }
    }

    public List<LineLink> getLinkedBreakends(){ return mLinkedBreakends; }
}
