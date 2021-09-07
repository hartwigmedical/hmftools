package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_BFB;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.types.ArmGroup;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class BfbFinder
{
    private CnDataLoader mCnDataLoader;

    public BfbFinder()
    {
        mCnDataLoader = null;
    }

    public void setCopyNumberAnalyser(CnDataLoader cnAnalyser) { mCnDataLoader = cnAnalyser; }

    public void analyseCluster(SvCluster cluster)
    {
        analyseCluster(cluster, false);
    }

    public void analyseCluster(SvCluster cluster, boolean reassess)
    {
        if(cluster.getFoldbacks().isEmpty() || cluster.getSvCount() == 1)
            return;

        final List<SvVarData> candidateDMSVs = Lists.newArrayList();

        double clusterMaxJcn = 0; // the max of the JCN min values

        double maxInfJcn = 0;

        for (SvVarData var : cluster.getSVs())
        {
            if (var.isSglBreakend())
                clusterMaxJcn = max(clusterMaxJcn, var.jcnMin() * 0.5); // in case the SGL is a disguised foldback
            else
                clusterMaxJcn = max(clusterMaxJcn, var.jcnMin());

            if (var.isInferredSgl())
                maxInfJcn = max(maxInfJcn, var.jcn());
        }

        // determine whether foldbacks and/or SGLs could explain the amplification, indicating a BFB process
        double sumFbJcn = 0;
        double maxFbJcn = 0;
        double foldbackCount = 0;
        double maxSglJcn = 0;

        for (final SvVarData var : cluster.getSVs())
        {
            if(candidateDMSVs.contains(var))
                continue;

            if(var.isFoldback())
            {
                sumFbJcn += var.isChainedFoldback() ? var.jcn() * 0.5 : var.jcn();
                foldbackCount += var.isChainedFoldback() ? 0.5 : 1;
                maxFbJcn = max(maxFbJcn, var.jcn());
            }
            else if(var.isSglBreakend())
            {
                // limit the SGLs to the top 2 by JCN
                maxSglJcn = max(var.jcn(), maxSglJcn);
            }
        }

        // determine the maximum potential JCN plausibly explained by BFB - taking the min of:
        // - 2x sum of FB ploidies + max INF /SGL JCN
        // - 6x max FB JCN
        // - HighestTelo/CentromereCN x 2^(FBCount + if(SGLPloidy > 10% * highest minPloidy,1,0))
        double maxArmEndCopyNumber = getMaxArmEndCopyNumber(cluster);
        double fbSglCount = foldbackCount + (maxSglJcn > 0.1 * clusterMaxJcn ? 1 : 0);
        double armBfbCopyNumber = maxArmEndCopyNumber * pow(2, fbSglCount);

        double maxBFBJcn = min(2 * sumFbJcn + maxSglJcn, min(6 * maxFbJcn, armBfbCopyNumber));

        if(maxBFBJcn > clusterMaxJcn && foldbackCount > 0)
        {
            LNX_LOGGER.debug(String.format("cluster(%s) BFB maxJCN(%.1f) plausible jcn(%.1f fb=%.1f sgl=%.1f arm=%.1f)",
                    cluster.id(), clusterMaxJcn, maxBFBJcn, sumFbJcn, maxSglJcn, armBfbCopyNumber));

            cluster.addAnnotation(CLUSTER_ANNOT_BFB);
        }
    }

    private double getMaxArmEndCopyNumber(final SvCluster cluster)
    {
        double maxCopyNumber = 0;

        for (ArmGroup armGroup : cluster.getArmGroups())
        {
            final TelomereCentromereCnData tcData = mCnDataLoader.getChrTeleCentroData().get(armGroup.chromosome());

            if(tcData != null)
            {
                double centromereCopyNumber = armGroup.arm() == P_ARM ? tcData.CentromerePArm : tcData.CentromereQArm;
                double telomereCopyNumber = armGroup.arm() == P_ARM ? tcData.TelomerePArm : tcData.TelomereQArm;
                maxCopyNumber = max(max(telomereCopyNumber, centromereCopyNumber), maxCopyNumber);
            }
        }

        return maxCopyNumber;
    }

}
