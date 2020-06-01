package com.hartwig.hmftools.linx.ext_compare;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvCluster;

import org.apache.commons.compress.utils.Lists;

public class CfChainClusterOverlap
{
    public final SvCluster Cluster;
    public final CfChain Chain;
    public final List<CfSvChainData> SharedSVs;

    private final String mHashId;

    public CfChainClusterOverlap(final SvCluster cluster, final CfChain chain)
    {
        Cluster = cluster;
        Chain = chain;
        SharedSVs = Lists.newArrayList();

        mHashId = formClusterChainId(Cluster.id(), Chain.ChainId);
    }

    public static String formClusterChainId(int clusterId, int chainId) { return String.format("%d_%d", clusterId, chainId); }

    public String hasdId() { return mHashId; }

}
