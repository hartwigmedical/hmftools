package com.hartwig.hmftools.linx.ext_compare;

import static com.hartwig.hmftools.linx.LinxOutput.SUBSET_SPLIT;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.linx.types.SvCluster;

import org.apache.commons.compress.utils.Lists;

public class CfChainClusterOverlap
{
    public final int Id;
    private final List<SvCluster> mClusters;
    private final List<CfChain> mChains;

    public CfChainClusterOverlap(int id)
    {
        Id = id;
        mClusters = Lists.newArrayList();
        mChains = Lists.newArrayList();
    }

    public final List<SvCluster> clusters() { return mClusters; }
    public final List<CfChain> chains() { return mChains; }

    public void addCluster(final SvCluster cluster)
    {
        // add in order
        for(int index = 0; index < mClusters.size(); ++index)
        {
            final SvCluster existingCluster = mClusters.get(index);

            if(existingCluster == cluster)
                return;

            if(cluster.getSvCount() > existingCluster.getSvCount())
            {
                mClusters.add(index, cluster);
                return;
            }
        }

        mClusters.add(cluster);
    }

    public void addChain(final CfChain chain)
    {
        for(int index = 0; index < mChains.size(); ++index)
        {
            final CfChain existingChain = mChains.get(index);

            if(existingChain == chain)
                return;

            if(chain.SVs.size() > existingChain.SVs.size())
            {
                mChains.add(index, chain);
                return;
            }
        }

        mChains.add(chain);
    }

    public static String header()
    {
        return "SampleId,OverlapGroupId,Clusters,Chains,TotalSVs,SharedSVs,ClusteredSvCount,SimpleSVs,SGLs,ChainedSvCount"
                + ",ClusterIds,ClusterCounts,ResolvedTypes,ChainIds,ChainCounts";
    }

    public String toCsv()
    {
        final StringJoiner csvOutput = new StringJoiner(",");

        csvOutput.add(String.valueOf(Id));
        csvOutput.add(String.valueOf(mClusters.size()));
        csvOutput.add(String.valueOf(mChains.size()));

        // count of SVs in both a chain and cluster in this group
        StringJoiner clusterIds = new StringJoiner(SUBSET_SPLIT);
        StringJoiner clusterCounts = new StringJoiner(SUBSET_SPLIT);
        StringJoiner resolvedTypes = new StringJoiner(SUBSET_SPLIT);

        int totalSVs = 0;
        int sglSVs = 0;
        int simpleSVs = 0;
        int clusteredSVs = 0;

        for(final SvCluster cluster : mClusters)
        {
            totalSVs += cluster.getSvCount();
            sglSVs += cluster.getSglBreakendCount();

            if(cluster.getSvCount() > 1)
            {
                clusteredSVs += cluster.getSvCount();
                clusterIds.add(String.valueOf(cluster.id()));
                clusterCounts.add(String.valueOf(cluster.getSvCount()));
                resolvedTypes.add(cluster.getResolvedType().toString());
            }
            else
            {
                ++simpleSVs;
            }
        }

        StringJoiner chainIds = new StringJoiner(SUBSET_SPLIT);
        StringJoiner chainCounts = new StringJoiner(SUBSET_SPLIT);

        int sharedSVs = 0;
        int chainedSVs = 0;
        for(final CfChain chain : mChains)
        {
            sharedSVs += chain.SVs.stream().filter(x -> x.getSvData() != null && x.getSvData().getCluster().getSvCount() > 1).count();
            chainedSVs += chain.SVs.size();

            chainIds.add(String.valueOf(chain.Id));
            chainCounts.add(String.valueOf(chain.SVs.size()));
        }

        csvOutput.add(String.valueOf(totalSVs));
        csvOutput.add(String.valueOf(sharedSVs));
        csvOutput.add(String.valueOf(clusteredSVs));
        csvOutput.add(String.valueOf(simpleSVs));
        csvOutput.add(String.valueOf(sglSVs));
        csvOutput.add(String.valueOf(chainedSVs));

        csvOutput.add(clusterIds.toString());
        csvOutput.add(clusterCounts.toString());
        csvOutput.add(resolvedTypes.toString());
        csvOutput.add(chainIds.toString());
        csvOutput.add(chainCounts.toString());

        return csvOutput.toString();
    }

    public String toString()
    {
        return String.format("%d: clusters(%d) chains(%d)", Id, mClusters.size(), mChains.size());
    }

}
