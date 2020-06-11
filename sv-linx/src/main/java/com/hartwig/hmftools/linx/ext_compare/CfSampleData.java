package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.linx.ext_compare.CfBreakendData.NO_ID;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class CfSampleData
{
    public final String SampleId;
    public final List<CfSvChainData> CfSvList;
    public final List<SvVarData> UnchainedSvList; // SVs not clustered by CF
    public final Map<Integer,CfChain> Chains;

    public final Map<String,CfChainClusterOverlap> LinkedClusterChains;

    public final Map<SvVarData,Integer> SvClusteringDistance;

    public CfSampleData(final String sampleId)
    {
        SampleId = sampleId;
        Chains = Maps.newHashMap();
        CfSvList = Lists.newArrayList();
        UnchainedSvList = Lists.newArrayList();
        LinkedClusterChains = Maps.newHashMap();
        SvClusteringDistance = Maps.newHashMap();
    }

    public void processNewSV(final CfSvChainData cfSvData)
    {
        // keep a cache of chains
        CfChain chain = Chains.get(cfSvData.ChainId);

        if(chain == null)
        {
            chain = new CfChain(cfSvData.ChainId);
            Chains.put(cfSvData.ChainId, chain);
        }

        chain.ChainSVs.add(cfSvData);

        if(cfSvData.svMatched())
        {
            UnchainedSvList.remove(cfSvData.getSvData());

            CfChainClusterOverlap clusterOverlap = LinkedClusterChains.get(cfSvData.clusterChainId());
            if(clusterOverlap == null)
            {
                clusterOverlap = new CfChainClusterOverlap(cfSvData.getSvData().getCluster(), chain);
                LinkedClusterChains.put(cfSvData.clusterChainId(), clusterOverlap);
            }

            clusterOverlap.SharedSVs.add(cfSvData);
        }
    }

    public void setDeletionBridgeData()
    {
        for(CfSvChainData cfSvData : CfSvList)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(cfSvData.DbBreakendIds[se] == NO_ID || cfSvData.getDbSvData()[se] != null) // non-existent or already matched
                    continue;

                boolean found = false;

                for(CfSvChainData otherSv : CfSvList)
                {
                    if(otherSv == cfSvData)
                        continue;

                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        if(cfSvData.DbBreakendIds[se] == otherSv.BreakpointIds[se2])
                        {
                            int dbLength = cfSvData.calcDbLength(otherSv, se, se2);
                            cfSvData.setDbSvData(se, otherSv.getSvData(), dbLength);
                            otherSv.setDbSvData(se2, cfSvData.getSvData(), dbLength);
                            found = true;
                            break;
                        }
                    }

                    if(found)
                        break;
                }
            }

        }
    }

    public void setSvClusterDistances(final SvCluster cluster)
    {
        if(cluster.getSvCount() == 1)
            return;

        for(int i = 0; i < cluster.getSvCount(); ++i)
        {
            final SvVarData var1 = cluster.getSV(i);

            final int otherSvId = Integer.parseInt(var1.getClusterReason().split(";", -1)[0].split("_")[1]);

            for(int j = 0; j < cluster.getSvCount(); ++j)
            {
                if(j == i)
                    continue;

                final SvVarData var2 = cluster.getSV(j);

                if(var2.id() == otherSvId)
                {
                    int distance = getProximity(var1, var2);

                    if(SvClusteringDistance.containsKey(var1))
                        SvClusteringDistance.put(var1, min(distance, SvClusteringDistance.get(var1)));
                    else
                        SvClusteringDistance.put(var1, distance);

                    if(SvClusteringDistance.containsKey(var2))
                        SvClusteringDistance.put(var2, min(distance, SvClusteringDistance.get(var2)));
                    else
                        SvClusteringDistance.put(var2, distance);

                    break;
                }
            }
        }
    }
}
