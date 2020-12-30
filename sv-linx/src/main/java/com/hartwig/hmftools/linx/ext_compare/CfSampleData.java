package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.linx.ext_compare.CfBreakendData.NO_ID;
import static com.hartwig.hmftools.linx.types.SvVarData.CR_DELIM;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class CfSampleData
{
    public final String SampleId;
    public final List<CfSvData> CfSvList;
    public final List<SvVarData> UnchainedSvList; // SVs not clustered by CF
    public final Map<Integer,CfChain> Chains;

    public final List<CfChainClusterOverlap> ClusterChainOverlaps;

    public final Map<SvVarData,Integer> SvProximityDistance;

    public CfSampleData(final String sampleId)
    {
        SampleId = sampleId;
        Chains = Maps.newHashMap();
        CfSvList = Lists.newArrayList();
        UnchainedSvList = Lists.newArrayList();
        ClusterChainOverlaps = Lists.newArrayList();
        SvProximityDistance = Maps.newHashMap();
    }

    public void processNewSV(final CfSvData cfSvData)
    {
        // keep a cache of chains
        CfChain chain = Chains.get(cfSvData.ChainId);

        if(chain == null)
        {
            chain = new CfChain(cfSvData.ChainId);
            Chains.put(cfSvData.ChainId, chain);
        }

        chain.addSV(cfSvData);

        if(cfSvData.svMatched())
        {
            UnchainedSvList.remove(cfSvData.getSvData());

            boolean clusterFound = false;

            for(CfChainClusterOverlap clusterOverlap : ClusterChainOverlaps)
            {
                boolean hasChain = clusterOverlap.chains().contains(chain);
                boolean hasCluster = clusterOverlap.clusters().contains(cfSvData.getSvData().getCluster());

                if(hasChain || hasCluster)
                {
                    clusterFound = true;

                    if(!hasCluster)
                        clusterOverlap.addCluster(cfSvData.getSvData().getCluster());

                    if(!hasChain)
                        clusterOverlap.addChain(chain);

                    break;
                }
            }

            if(!clusterFound)
            {
                CfChainClusterOverlap clusterOverlap = new CfChainClusterOverlap(ClusterChainOverlaps.size());
                clusterOverlap.addCluster(cfSvData.getSvData().getCluster());
                clusterOverlap.addChain(chain);
                ClusterChainOverlaps.add(clusterOverlap);
            }
        }

        reconcileOverlaps();
    }

    private void reconcileOverlaps()
    {
        if(ClusterChainOverlaps.size() == 1)
            return;

        int index = 0;
        while(index < ClusterChainOverlaps.size())
        {
            CfChainClusterOverlap clusterOverlap = ClusterChainOverlaps.get(index);

            boolean linkFound = false;

            for(int i = index + 1; i < ClusterChainOverlaps.size(); ++i)
            {
                CfChainClusterOverlap otherOverlap = ClusterChainOverlaps.get(i);

                if(clusterOverlap.clusters().stream().anyMatch(x -> otherOverlap.clusters().contains(x))
                || clusterOverlap.chains().stream().anyMatch(x -> otherOverlap.chains().contains(x)))
                {
                    linkFound = true;

                    otherOverlap.clusters().forEach(x -> clusterOverlap.addCluster(x));
                    otherOverlap.chains().forEach(x -> clusterOverlap.addChain(x));
                    ClusterChainOverlaps.remove(i);
                    break;
                }
            }

            if(!linkFound)
            {
                ++index;
            }
        }
    }

    public final CfChainClusterOverlap findClusterChainOverlap(final CfChain chain)
    {
        return ClusterChainOverlaps.stream().filter(x -> x.chains().contains(chain)).findFirst().orElse(null);
    }

    public void setDeletionBridgeData()
    {
        for(CfSvData cfSvData : CfSvList)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(cfSvData.DbBreakendIds[se] == NO_ID || cfSvData.getDbSvData()[se] != null) // non-existent or already matched
                    continue;

                boolean found = false;

                for(CfSvData otherSv : CfSvList)
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

    public void setChainLinks()
    {
        Chains.values().forEach(x -> x.buildLinks());
    }

    public void setSvClusterDistances(final SvCluster cluster)
    {
        // set the distance to the first SV marked in the clustering reasons
        if(cluster.getSvCount() == 1)
            return;

        for(int i = 0; i < cluster.getSvCount(); ++i)
        {
            final SvVarData var1 = cluster.getSV(i);

            final int otherSvId = Integer.parseInt(var1.getClusterReason().split(ITEM_DELIM, -1)[0].split(CR_DELIM)[1]);

            for(int j = 0; j < cluster.getSVs().size(); ++j)
            {
                if(j == i)
                    continue;

                final SvVarData var2 = cluster.getSV(j);

                if(var2.id() == otherSvId)
                {
                    int distance = getProximity(var1, var2);

                    if(SvProximityDistance.containsKey(var1))
                        SvProximityDistance.put(var1, min(distance, SvProximityDistance.get(var1)));
                    else
                        SvProximityDistance.put(var1, distance);

                    if(SvProximityDistance.containsKey(var2))
                        SvProximityDistance.put(var2, min(distance, SvProximityDistance.get(var2)));
                    else
                        SvProximityDistance.put(var2, distance);

                    break;
                }
            }
        }
    }

    public void setUnclusteredDistances(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        // sets distance to nearest SVs for unclustered SVs
        for(CfSvData cfSvData : CfSvList)
        {
            final SvVarData var = cfSvData.getSvData();

            if(var == null || var.isSglBreakend())
                continue;

            if(var.getCluster().getSvCount() > 1)
                continue;

            // find the nearest SV distance for this unclustered SV
            int minDistance = -1;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = var.getBreakend(se);

                final List<SvBreakend> breakendList = chrBreakendMap.get(breakend.chromosome());

                for(final SvBreakend otherBreakend : breakendList)
                {
                    if(otherBreakend.getSV() == breakend.getSV())
                        continue;

                    int distance = abs(breakend.position() - otherBreakend.position());

                    if(minDistance == -1 || distance < minDistance)
                    {
                        minDistance = distance;
                    }
                }
            }

            if(minDistance >= 0)
                SvProximityDistance.put(var, minDistance);
        }
    }
}
