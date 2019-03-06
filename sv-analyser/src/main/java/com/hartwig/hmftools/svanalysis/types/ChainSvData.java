package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;

public class ChainSvData
{
    public final SvVarData SV;
    public SvVarData DualFoldbackOtherSV;

    public int ReplicationCount;

    public List<SvVarData> OverlappedFoldbacks;
    public List<SvBreakend> BreakendsFacingStart;
    public List<SvBreakend> BreakendsFacingEnd;

    public String Type;

    public static final String CHAIN_SV_TYPE_FOLDBACK = "Foldback";
    public static final String CHAIN_SV_TYPE_COMPLEX_DUP = "ComplexDUP";
    public static final String CHAIN_SV_TYPE_COMPLEX_INV = "ComplexINV";
    public static final String CHAIN_SV_TYPE_SIMPLE = "Simple";

    public ChainSvData(final SvVarData sv)
    {
        SV = sv;
        DualFoldbackOtherSV= null;
        OverlappedFoldbacks = Lists.newArrayList();
        BreakendsFacingStart = Lists.newArrayList();
        BreakendsFacingEnd = Lists.newArrayList();
        ReplicationCount = max(sv.getReplicatedCount(), 1);
        Type = "";
    }

    public List<SvBreakend> getFacingBreakends(boolean isStart)
    {
        return isStart ? BreakendsFacingStart : BreakendsFacingEnd;
    }

    public final SvBreakend getBreakend(boolean useStart)
    {
        if(Type != CHAIN_SV_TYPE_FOLDBACK || DualFoldbackOtherSV == null)
        {
            return SV.getBreakend(useStart);
        }

        SvBreakend thisBreakend;
        SvBreakend otherBreakend;

        if(SV.getFoldbackBreakend(true) != null)
        {
            thisBreakend = SV.getBreakend(true);
            otherBreakend = SV.getFoldbackBreakend(true);
        }
        else
        {
            thisBreakend = SV.getBreakend(false);
            otherBreakend = SV.getFoldbackBreakend(false);
        }

        return (useStart) == (thisBreakend.position() < otherBreakend.position()) ? thisBreakend : otherBreakend;

    }

}
