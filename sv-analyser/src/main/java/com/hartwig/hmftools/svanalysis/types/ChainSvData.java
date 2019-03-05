package com.hartwig.hmftools.svanalysis.types;

import static java.lang.Math.max;

import java.util.List;

import com.google.common.collect.Lists;

public class ChainSvData
{
    public final SvVarData SV;
    public int ReplicationCount;

    public List<SvVarData> OverlappedFoldbacks;
    public List<SvBreakend> BreakendsFacedStart;
    public List<SvBreakend> BreakendsFacedEnd;

    public String Type;

    public static final String CHAIN_SV_TYPE_FOLDBACK = "Foldback";
    public static final String CHAIN_SV_TYPE_COMPLEX_DUP = "ComplexDUP";
    public static final String CHAIN_SV_TYPE_COMPLEX_INV = "ComplexINV";
    public static final String CHAIN_SV_TYPE_SIMPLE = "Simple";

    public ChainSvData(final SvVarData sv)
    {
        SV = sv;
        OverlappedFoldbacks = Lists.newArrayList();
        BreakendsFacedStart = Lists.newArrayList();
        BreakendsFacedEnd = Lists.newArrayList();
        ReplicationCount = max(sv.getReplicatedCount(), 1);
        Type = "";
    }

}
