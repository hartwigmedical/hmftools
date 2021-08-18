package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.linx.drivers.DriverAmpData;
import com.hartwig.hmftools.linx.drivers.OpposingSegment;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;

public class AmpliconClusterMatch
{
    public static final double MIN_AMP_PERCENT_VS_MAX = 0.33;

    public static List<DriverAmpData> findAmplifyingClusters(
            final ChrBaseRegion ampRegion, final List<SvBreakend> breakendList)
    {
        if (breakendList == null || breakendList.isEmpty())
            return Lists.newArrayList();

        final Map<SvCluster,DriverAmpData> clusterAmpData = Maps.newHashMap();

        // for each cluster, calculate the location and highest CN within the amplified region
        // then calculate the net contribution towards this point but no further

        // double ampCopyNumber = calcRegionCopyNumber(ampRegion, breakendList);

        // sum up breakend ploidies from telomere to centromere for the gene in question net off ploidy within a cluster
        int startIndex = 0;
        int endIndex = breakendList.size() - 1;

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
                if ((traverseUp && breakend.position() > ampRegion.end()) || (!traverseUp && breakend.position() < ampRegion.start()))
                    break;

                final SvCluster cluster = breakend.getCluster();

                if (cluster.hasLinkingLineElements() || processedClusters.contains(cluster))
                    continue;

                if(cluster.getSvCount() == 1 && cluster.getSV(0).type() == DEL) // explicitly disallow
                    continue;

                processedClusters.add(cluster);

                // proceed from this point until the start of the gene
                DriverAmpData ampData = checkClusterForAmplification(ampRegion, breakendList, breakend, traverseUp, opposingSegments);

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
            return Lists.newArrayList();

        return clusterAmpData.values().stream().collect(Collectors.toList());
    }

    private static DriverAmpData checkClusterForAmplification(
            final ChrBaseRegion ampRegion, final List<SvBreakend> breakendList, final SvBreakend startBreakend,
            boolean traverseUp, final List<OpposingSegment> opposingSegments)
    {
        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);

        final SvCluster targetCluster = startBreakend.getCluster();

        boolean inSegment = true;
        double segmentStartCopyNumber = startCopyNumber;
        double segmentMaxCopyNumber = startBreakend.copyNumber();
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

            if ((traverseUp && breakend.position() > ampRegion.end()) || (!traverseUp && breakend.position() < ampRegion.start()))
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

                    inSegment = false;
                }

                OpposingSegment opposingSegment = findOrCreateOpposingSegment(
                        opposingSegments, breakendList, breakend, traverseUp, traverseUp ? ampRegion.start() : ampRegion.end());

                if(netClusterCNChange > 0 && opposingSegment != null && opposingSegment.remainingCNChange() > 0)
                {
                    if(opposingSegment.remainingCNChange() > netClusterCNChange)
                    {
                        LNX_LOGGER.trace("ampRegion({}) cluster({}) netCN({}) cancelled by breakend({}) cnChange(orig={} remain={})",
                                ampRegion, targetCluster.id(), formatJcn(netClusterCNChange),
                                breakend, formatJcn(breakend.copyNumberChange()), formatJcn(opposingSegment.remainingCNChange()));

                        opposingSegment.reduceCNChange(netClusterCNChange);
                        netClusterCNChange = 0;
                    }
                    else
                    {
                        LNX_LOGGER.trace("ampRegion({}) cluster({}) netCN({}) reducing by breakend({}) cnChange({})",
                                ampRegion, targetCluster.id(), formatJcn(netClusterCNChange),
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
                    segmentMaxCopyNumber = breakend.copyNumber();
                    segStartBreakend = breakend;
                    inSegment = true;
                }
                else
                {
                    segmentMaxCopyNumber = max(segmentMaxCopyNumber, breakend.copyNumber());
                }
            }
        }

        double clusterCNChange = netClusterCNChange;

        if(inSegment)
        {
            if(segmentMaxCopyNumber > segmentStartCopyNumber)
                clusterCNChange += segmentMaxCopyNumber - segmentStartCopyNumber;

            LNX_LOGGER.trace("ampRegion({}) cluster({}) open segment startCN({}) net({}) start breakend({})",
                    ampRegion, targetCluster.id(), formatJcn(segmentStartCopyNumber), formatJcn(clusterCNChange), segStartBreakend);
        }

        if(clusterCNChange > 0 && !copyNumbersEqual(clusterCNChange, 0) && !copyNumbersEqual(startCopyNumber + clusterCNChange, startCopyNumber))
        {
            LNX_LOGGER.debug("ampRegion({}) cluster({}) copy number gain({}) vs startCN({}) segments({}: CN={}) traversal({})",
                    ampRegion, targetCluster.id(), formatJcn(clusterCNChange), formatJcn(startCopyNumber),
                    segmentCount, formatJcn(netClusterCNChange), traverseUp ? "up" : "down");

            return new DriverAmpData(targetCluster, traverseUp, breakendCount, segmentCount, startCopyNumber, clusterCNChange);
        }

        return null;
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

        opposingSegments.add(opposingSegment);
        return opposingSegment;
    }

    private static double calcRegionCopyNumber(final ChrBaseRegion region, final List<SvBreakend> breakendList)
    {
        double regionCopyNumber = -1;

        for(int i = 0; i < breakendList.size(); ++i)
        {
            SvBreakend breakend = breakendList.get(i);
            SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

            if(breakend.position() < region.start())
            {
                if(nextBreakend == null)
                {
                    return breakend.copyNumber();
                }
                else
                {
                    if(nextBreakend.position() < region.start())
                        continue;

                    // next breakend is past the start of the target region
                    regionCopyNumber = breakend.getCopyNumber(false);
                }
            }
            else if(breakend.position() > region.end())
            {
                break;
            }
            else
            {
                regionCopyNumber = max(breakend.copyNumber(), regionCopyNumber);
            }
        }

        return regionCopyNumber;
    }
}
