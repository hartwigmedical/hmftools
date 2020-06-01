package com.hartwig.hmftools.linx.ext_compare;

import java.util.List;

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

}
