package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.isShortArmChromosome;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.common.linx.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;

public class AmplificationDrivers
{
    private final DriverDataCache mDataCache;

    public AmplificationDrivers(final DriverDataCache dataCache)
    {
        mDataCache = dataCache;
    }

    public void annotateAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        LNX_LOGGER.debug("gene({}) chromosome({}) position({} -> {}) minCN({})",
                dgData.GeneInfo.GeneName, dgData.GeneInfo.Chromosome, dgData.TransData.TransStart, dgData.TransData.TransEnd,
                formatJcn(dgData.CopyNumberRegion.MinCopyNumber));

        checkChromosomeAmplification(dgData);
        checkClusterAmplification(dgData, breakendList);

        if(dgData.getEvents().isEmpty())
        {
            LNX_LOGGER.debug("gene({}) AMP gain no event found", dgData.GeneInfo.GeneName);

            // otherwise no event
            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            dgData.addEvent(event);
        }
    }

    private void checkChromosomeAmplification(final DriverGeneData dgData)
    {
        // check for an arm or whole chromosome amplified above the sample ploidy
        final TelomereCentromereCnData tcData = mDataCache.CopyNumberData.getChrTeleCentroData().get(dgData.GeneInfo.Chromosome);

        if(tcData == null)
            return;

        double telomereCopyNumber = min(tcData.TelomerePArm, tcData.TelomereQArm);

        double centromereCopyNumber = min(tcData.CentromerePArm, tcData.CentromereQArm);

        double chromosomeCopyNumber = min(centromereCopyNumber, telomereCopyNumber);

        double samplePloidy = mDataCache.samplePloidy();

        if(!isShortArmChromosome(dgData.GeneInfo.Chromosome))
        {
            if (chromosomeCopyNumber / samplePloidy > 2)
            {
                LNX_LOGGER.debug("gene({}) AMP gain from chromosome chChange({} telo={} centro={}) vs samplePloidy({})",
                        dgData.GeneInfo.GeneName, formatJcn(chromosomeCopyNumber),
                        formatJcn(telomereCopyNumber), formatJcn(centromereCopyNumber), formatPloidy(samplePloidy));

                DriverGeneEvent event = new DriverGeneEvent(GAIN_CHR);
                event.setCopyNumberGain(chromosomeCopyNumber - samplePloidy);
                dgData.addEvent(event);
            }
        }

        double centromereCNChange = dgData.Arm == P_ARM ?
                min(tcData.CentromerePArm, tcData.TelomerePArm) - tcData.CentromereQArm
                : min(tcData.CentromereQArm, tcData.TelomereQArm) - tcData.CentromerePArm;

        if (centromereCNChange > 0 && !copyNumbersEqual(centromereCNChange, 0))
        {
            LNX_LOGGER.debug("gene({}) AMP gain from arm({}) cnChange across centromere({} -> {} = {})",
                    dgData.GeneInfo.GeneName, dgData.Arm, formatJcn(tcData.CentromerePArm), formatJcn(tcData.CentromereQArm),
                    formatJcn(centromereCNChange));

            DriverGeneEvent event = new DriverGeneEvent(GAIN_ARM);
            event.setCopyNumberGain(centromereCNChange);
            dgData.addEvent(event);
        }
    }

    private static OpposingSegment findOrCreateOpposingSegment(
            List<OpposingSegment> opposingSegments, final List<SvBreakend> breakendList, final SvBreakend startBreakend,
            boolean traverseUp, int stopPosition)
    {
        OpposingSegment opposingSegment = opposingSegments.stream()
                .filter(x -> x.Breakends.contains(startBreakend))
                .findFirst().orElse(null);

        if(opposingSegment != null)
            return opposingSegment;

        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);
        double endCopyNumber = startBreakend.getCopyNumber(!traverseUp);

        int index = startBreakend.getChrPosIndex();

        List<SvBreakend> segmentBreakends = Lists.newArrayList(startBreakend);

        while (true)
        {
            index += traverseUp ? 1 : -1;

            if (index < 0 || index >= breakendList.size())
                break;

            SvBreakend breakend = breakendList.get(index);

            if ((traverseUp && breakend.position() > stopPosition) || (!traverseUp && breakend.position() < stopPosition))
            {
                break;
            }

            if (breakend.getCluster() == startBreakend.getCluster())
            {
                segmentBreakends.add(breakend);
                endCopyNumber = breakend.getCopyNumber(!traverseUp);
            }
            else
            {
                break;
            }
        }

        double netClusterCNChange = endCopyNumber - startCopyNumber;

        if(netClusterCNChange > 0 || copyNumbersEqual(netClusterCNChange, 0))
        {
            opposingSegment = new OpposingSegment(startBreakend.getCluster(), segmentBreakends, 0);
            opposingSegments.add(opposingSegment);
            return null;
        }

        opposingSegment = new OpposingSegment(startBreakend.getCluster(), segmentBreakends, -netClusterCNChange);

        LNX_LOGGER.trace("added opposing segment: {}", opposingSegment);

        opposingSegments.add(opposingSegment);
        return opposingSegment;
    }

    private static DriverAmpData checkClusterForAmplification(
            final DriverGeneData dgData, final List<SvBreakend> breakendList, final SvBreakend startBreakend, boolean traverseUp,
            List<OpposingSegment> opposingSegments)
    {
        int transStart = dgData.TransData.TransStart;
        int transEnd = dgData.TransData.TransEnd;

        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);

        final SvCluster targetCluster = startBreakend.getCluster();

        boolean inSegment = true;
        double segmentStartCopyNumber = startCopyNumber;
        double netClusterCNChange = 0;
        int segmentCount = 0;
        int breakendCount = 0;

        int index = startBreakend.getChrPosIndex();
        SvBreakend breakend = null;
        SvBreakend segStartBreakend = startBreakend;

        while (true)
        {
            if (breakend != null)
                index += traverseUp ? 1 : -1;

            if(index < 0 || index >= breakendList.size())
                break;

            breakend = breakendList.get(index);

            if ((traverseUp && breakend.position() > transStart) || (!traverseUp && breakend.position() < transEnd))
                break;

            final SvCluster cluster = breakend.getCluster();

            if(cluster != targetCluster)
            {
                if(inSegment)
                {
                    // record details of this segment
                    double endCopyNumber = 0;
                    if(breakend.arm().equals(segStartBreakend.arm()))
                    {
                        endCopyNumber = breakend.getCopyNumber(traverseUp);
                    }
                    else
                    {
                        // take the end copy number from the centromere
                        int prevIndex = index + (traverseUp ? -1 : 1);
                        final SvBreakend prevBreakend = breakendList.get(prevIndex);
                        endCopyNumber = prevBreakend.getCopyNumber(!traverseUp);
                    }

                    double segmentCNChange = endCopyNumber - segmentStartCopyNumber;
                    netClusterCNChange += segmentCNChange;

                    LNX_LOGGER.trace("gene({}) cluster({}) adding segment CN({} -> {} chg={}) net({}) breakends({} -> {})",
                            dgData.GeneInfo.GeneName, targetCluster.id(), formatJcn(segmentStartCopyNumber), formatJcn(endCopyNumber),
                            formatJcn(segmentCNChange), formatJcn(netClusterCNChange), segStartBreakend, breakend);

                    inSegment = false;
                }

                OpposingSegment opposingSegment = findOrCreateOpposingSegment(
                        opposingSegments, breakendList, breakend, traverseUp, traverseUp ? transStart : transEnd);

                if(netClusterCNChange > 0 && opposingSegment != null && opposingSegment.remainingCNChange() > 0)
                {
                    if(opposingSegment.remainingCNChange() > netClusterCNChange)
                    {
                        LNX_LOGGER.trace("gene({}) cluster({}) netCN({}) cancelled by breakend({}) cnChange(orig={} remain={})",
                                dgData.GeneInfo.GeneName, targetCluster.id(), formatJcn(netClusterCNChange),
                                breakend, formatJcn(breakend.copyNumberChange()), formatJcn(opposingSegment.remainingCNChange()));

                        opposingSegment.reduceCNChange(netClusterCNChange);
                        netClusterCNChange = 0;
                    }
                    else
                    {
                        LNX_LOGGER.trace("gene({}) cluster({}) netCN({}) reducing by breakend({}) cnChange({})",
                                dgData.GeneInfo.GeneName, targetCluster.id(), formatJcn(netClusterCNChange),
                                breakend, formatJcn(breakend.copyNumberChange()));

                        netClusterCNChange -= opposingSegment.remainingCNChange();
                        opposingSegment.zeroCNChange(); // keep so it won't be registered again but zero out
                    }
                }
            }
            else if(cluster == targetCluster)
            {
                ++breakendCount;

                if(!inSegment)
                {
                    segmentStartCopyNumber = breakend.getCopyNumber(traverseUp);
                    segStartBreakend = breakend;
                    inSegment = true;
                }
            }
        }

        double geneMinCopyNumber = dgData.CopyNumberRegion.MinCopyNumber;

        double clusterCNChange = netClusterCNChange;

        if(inSegment)
        {
            clusterCNChange += geneMinCopyNumber - segmentStartCopyNumber;

            LNX_LOGGER.trace("gene({}) cluster({}) open segment startCN({}) net({}) start breakend({})",
                    dgData.GeneInfo.GeneName, targetCluster.id(), formatJcn(segmentStartCopyNumber), formatJcn(clusterCNChange),
                    segStartBreakend);
        }

        if(clusterCNChange > 0 && !copyNumbersEqual(clusterCNChange, 0) && !copyNumbersEqual(startCopyNumber + clusterCNChange, startCopyNumber))
        {
            LNX_LOGGER.debug("gene({}) cluster({}) copy number gain({}) vs startCN({}) segments({}: CN={}) traversal({})",
                    dgData.GeneInfo.GeneName, targetCluster.id(), formatJcn(clusterCNChange), formatJcn(startCopyNumber),
                    segmentCount, formatJcn(netClusterCNChange), traverseUp ? "up" : "down");

            return new DriverAmpData(targetCluster, traverseUp, breakendCount, segmentCount, startCopyNumber, clusterCNChange);
        }

        return null;
    }

    private static final double MIN_AMP_PERCENT_VS_MAX = 0.33;

    private static void checkClusterAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        if (breakendList == null || breakendList.isEmpty())
            return;

        // sum up breakend ploidies from telomere to centromere for the gene in question net off ploidy within a cluster
        int transStart = dgData.TransData.TransStart;
        int transEnd = dgData.TransData.TransEnd;

        int startIndex = 0;
        int endIndex = breakendList.size() - 1;

        Map<SvCluster,DriverAmpData> clusterAmpData = Maps.newHashMap();

        for(int i = 0; i <= 1; ++i)
        {
            boolean traverseUp = (i == 0);

            List<SvCluster> processedClusters = Lists.newArrayList();
            List<OpposingSegment> opposingSegments = Lists.newArrayList();

            int index = traverseUp ? startIndex : endIndex;
            SvBreakend breakend = null;

            while (true)
            {
                if (breakend != null)
                    index += traverseUp ? 1 : -1;

                if (index < startIndex || index > endIndex)
                    break;

                breakend = breakendList.get(index);

                // make note of any breakend preceding the gene in the direction of traversal
                if ((traverseUp && breakend.position() > transEnd) || (!traverseUp && breakend.position() < transStart))
                    break;

                final SvCluster cluster = breakend.getCluster();

                if (cluster.hasLinkingLineElements() || processedClusters.contains(cluster))
                    continue;

                if(cluster.getSvCount() == 1 && cluster.getSV(0).type() == DEL) // explicitly disallow
                    continue;

                processedClusters.add(cluster);

                // proceed from this point until the start of the gene
                DriverAmpData ampData = checkClusterForAmplification(dgData, breakendList, breakend, traverseUp, opposingSegments);

                if(ampData == null)
                    continue;

                DriverAmpData existingAmpData = clusterAmpData.get(cluster);

                if(existingAmpData == null || existingAmpData.NetCNChange < ampData.NetCNChange)
                {
                    clusterAmpData.put(cluster, ampData);
                }
            }
        }

        if(clusterAmpData.isEmpty())
            return;

        // from the identified AMP clusters, find the max and percentage contributions of each
        double maxCnChange = clusterAmpData.values().stream().mapToDouble(x -> x.NetCNChange).max().getAsDouble();

        int index = 0;
        while(index < dgData.getEvents().size())
        {
            DriverGeneEvent event = dgData.getEvents().get(index);
            double armChrGain = event.getCopyNumberGain();

            if(armChrGain < MIN_AMP_PERCENT_VS_MAX * maxCnChange)
            {
                dgData.getEvents().remove(index);
            }
            else
            {
                maxCnChange = max(maxCnChange, armChrGain);
                ++index;
            }
        }

        for(Map.Entry<SvCluster,DriverAmpData> entry : clusterAmpData.entrySet())
        {
            final SvCluster cluster = entry.getKey();
            final DriverAmpData ampData = entry.getValue();

            // take any at least 20% of the largest contributing cluster
            if(ampData.NetCNChange / maxCnChange < MIN_AMP_PERCENT_VS_MAX)
                continue;

            LNX_LOGGER.debug("gene({}) cluster({}) adding AMP data: {}",
                    dgData.GeneInfo.GeneName, cluster.id(), ampData);

            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            event.setCluster(cluster);
            event.setCopyNumberGain(ampData.NetCNChange);
            event.setAmpData(ampData);

            event.setSvInfo(String.format("%s;%d;%d;%.1f",
                    ampData.TraverseUp, ampData.BreakendCount, ampData.SegmentCount, ampData.StartCopyNumber));
            dgData.addEvent(event);
        }
    }
}
