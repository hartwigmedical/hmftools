package com.hartwig.hmftools.linx.types;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.chaining.SvChain;

public class DoubleMinuteData
{
    public final SvCluster Cluster;
    public final List<SvVarData> SVs;

    public double MaxBFBJcn;
    public double MinAdjacentMARatio;

    public final List<SvVarData> CandidateSVs;
    public final List<SvChain> Chains;
    public boolean FullyChained;

    public DoubleMinuteData(final SvCluster cluster, final List<SvVarData> svList)
    {
        Cluster = cluster;
        SVs = svList;

        MinAdjacentMARatio = 0;
        MaxBFBJcn = 0;

        Chains = Lists.newArrayList();
        CandidateSVs = Lists.newArrayList();
        FullyChained = false;
    }

}
