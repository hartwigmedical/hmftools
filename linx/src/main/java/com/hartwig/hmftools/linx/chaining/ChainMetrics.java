package com.hartwig.hmftools.linx.chaining;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_EXTERNAL;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_INTERNAL;
import static com.hartwig.hmftools.linx.types.LinkedPair.LOCATION_TYPE_REMOTE;

import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.LinkedPair;

public class ChainMetrics
{
    public int InternalTIs;
    public int InternalTICnGain;
    public int InternalShortTIs;
    public int ExternalTIs;
    public int ExternalTICnGain;
    public int ExternalShortTIs;
    public int OverlappingTIs;

    public int ChainEndsFace;
    public int ChainEndsAway;

    public ChainMetrics()
    {
        InternalTIs = 0;
        InternalTICnGain = 0;
        InternalShortTIs = 0;
        ExternalTIs = 0;
        ExternalTICnGain = 0;
        ExternalShortTIs = 0;
        OverlappingTIs = 0;
        ChainEndsFace = 0;
        ChainEndsAway = 0;
    }

    public void add(final ChainMetrics other)
    {
        InternalTIs += other.InternalTIs;
        InternalTICnGain += other.InternalTICnGain;
        InternalShortTIs += other.InternalShortTIs;
        ExternalTIs += other.ExternalTIs;
        ExternalTICnGain += other.ExternalTICnGain;
        ExternalShortTIs += other.ExternalShortTIs;
        OverlappingTIs += other.OverlappingTIs;
        ChainEndsFace += other.ChainEndsFace;
        ChainEndsAway += other.ChainEndsAway;
    }

    public static ChainMetrics extractChainMetrics(final SvChain chain)
    {
        final SvBreakend chainStart = chain.getOpenBreakend(true);
        final SvBreakend chainEnd = chain.getOpenBreakend(false);

        ChainMetrics metrics = new ChainMetrics();

        if(chainStart == null || chainEnd == null)
            return metrics;

        final SvBreakend lowerBreakend = chainStart.position() < chainEnd.position() ? chainStart : chainEnd;
        final SvBreakend upperBreakend = chainStart == lowerBreakend ? chainEnd : chainStart;

        boolean startEndSameArm = chainEnd.getChrArm().equals(chainStart.getChrArm());

        if(startEndSameArm)
        {
            if(lowerBreakend.orientation() == 1 && upperBreakend.orientation() == -1)
            {
                ++metrics.ChainEndsAway;
            }
            else if(lowerBreakend.orientation() == -1 && upperBreakend.orientation() == 1)
            {
                ++metrics.ChainEndsFace;
            }
        }

        for(final LinkedPair pair : chain.getLinkedPairs())
        {
            if(pair.first().type() == SGL || pair.second().type() == SGL)
                continue;

            if(pair.locationType() == LOCATION_TYPE_INTERNAL)
            {
                ++metrics.InternalTIs;

                if(pair.baseLength() <= SHORT_TI_LENGTH)
                    ++metrics.InternalShortTIs;

                if(pair.hasCopyNumberGain())
                    ++metrics.InternalTICnGain;
            }
            else if(pair.locationType() == LOCATION_TYPE_REMOTE || pair.locationType() == LOCATION_TYPE_EXTERNAL)
            {
                ++metrics.ExternalTIs;

                if(pair.baseLength() <= SHORT_TI_LENGTH)
                    ++metrics.ExternalShortTIs;

                if(pair.hasCopyNumberGain())
                    ++metrics.ExternalTICnGain;
            }

            if(pair.overlapCount() > 0)
                ++metrics.OverlappingTIs;
        }

        return metrics;
    }

}
