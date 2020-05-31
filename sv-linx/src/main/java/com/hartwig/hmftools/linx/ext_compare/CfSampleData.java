package com.hartwig.hmftools.linx.ext_compare;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class CfSampleData
{
    public final String SampleId;
    public final List<CfSvChainData> CfSvList;
    public final List<SvVarData> UnchainedSvList;
    public final Map<Integer,CfChain> Chains;

    public CfSampleData(final String sampleId)
    {
        SampleId = sampleId;
        Chains = Maps.newHashMap();
        CfSvList = Lists.newArrayList();
        UnchainedSvList = Lists.newArrayList();
    }
}
