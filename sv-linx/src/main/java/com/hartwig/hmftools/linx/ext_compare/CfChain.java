package com.hartwig.hmftools.linx.ext_compare;

import java.util.List;

import com.hartwig.hmftools.linx.types.SvCluster;

import org.apache.commons.compress.utils.Lists;

public class CfChain
{
    public final int ChainId;
    public final List<CfSvChainData> ChainSVs;

    public CfChain(int chainId)
    {
        ChainId = chainId;
        ChainSVs = Lists.newArrayList();
    }

    public int getSharedSvCount(final CfSvChainData cfSvData)
    {
        final SvCluster svCluster = cfSvData.getSvData().getCluster();
        return (int)ChainSVs.stream().filter(x -> x.getSvData().getCluster() == svCluster).count();
    }

}
