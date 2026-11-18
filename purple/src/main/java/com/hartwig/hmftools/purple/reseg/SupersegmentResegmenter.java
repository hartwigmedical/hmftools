package com.hartwig.hmftools.purple.reseg;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static java.lang.System.arraycopy;

import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_BAF_WEIGHT_THRESHOLD;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_SUPERSEGMENT_SIZE_FOR_BRUTE_FORCE;
import static com.hartwig.hmftools.purple.region.ObservedRegion.fromOther;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class SupersegmentResegmenter
{
    public static List<ObservedRegion> resegment(final List<ObservedRegion> members, double segmentationPenalty)
    {
        // calculate optimal segmentation based on the calculated tumor ratio penalty
        if(members.size() > RESEG_MAX_SUPERSEGMENT_SIZE_FOR_BRUTE_FORCE)
            return splitLargeSupersegment(members, segmentationPenalty);

        List<List<ObservedRegion>> bestPartition = findBestPartition(members, segmentationPenalty);

        if(bestPartition == null)
            return List.of(aggregateRegions(members));

        List<ObservedRegion> result = Lists.newArrayList();

        for(List<ObservedRegion> subsegment : bestPartition)
        {
            result.add(aggregateRegions(subsegment));
        }

        return result;
    }

    private static List<ObservedRegion> splitLargeSupersegment(final List<ObservedRegion> members, double segmentationPenalty)
    {
        double maxRatio = members.stream().mapToDouble(m -> m.observedTumorRatio()).max().orElse(0);
        double minRatio = members.stream().mapToDouble(m -> m.observedTumorRatio()).min().orElse(0);

        if(maxRatio - minRatio > sqrt(segmentationPenalty / 2))
        {
            int cutoffIndex = findLargestAdjacentDiffIndex(members);

            List<ObservedRegion> left = new ArrayList<>(members.subList(0, cutoffIndex));
            List<ObservedRegion> right = new ArrayList<>(members.subList(cutoffIndex, members.size()));

            List<ObservedRegion> result = new ArrayList<>();
            result.addAll(resegment(left, segmentationPenalty));
            result.addAll(resegment(right, segmentationPenalty));
            return result;
        }

        return List.of(aggregateRegions(members));
    }

    private static int findLargestAdjacentDiffIndex(final List<ObservedRegion> members)
    {
        int bestIndex = 1;
        double bestDiff = -1;

        for(int i = 1; i < members.size(); i++)
        {
            double diff = abs(members.get(i).observedTumorRatio() - members.get(i - 1).observedTumorRatio());

            if(diff >= bestDiff)
            {
                bestDiff = diff;
                bestIndex = i;
            }
        }

        return bestIndex;
    }

    private static List<List<ObservedRegion>> findBestPartition(final List<ObservedRegion> members, double segmentationPenalty)
    {
        // returns null if no partition beats the unpartitioned (0-cut) baseline - caller should aggregate members as one
        int n = members.size();

        double lowestPenalty = deviationPenalty(members);
        int[] bestPartitionPoints = null;

        for(int k = 1; k < n; k++)
        {
            double partitionCost = k * segmentationPenalty;

            if(partitionCost >= lowestPenalty)
                break;

            for(int[] cuts : combinationsOf(n - 1, k))
            {
                int[] padded = pad(cuts, n);
                double deviationCost = 0;

                for(int j = 0; j < padded.length - 1; j++)
                {
                    List<ObservedRegion> subDf = members.subList(padded[j], padded[j + 1]);
                    deviationCost += deviationPenalty(subDf);

                    if(partitionCost + deviationCost >= lowestPenalty)
                        break;
                }

                if(partitionCost + deviationCost < lowestPenalty)
                {
                    bestPartitionPoints = padded;
                    lowestPenalty = partitionCost + deviationCost;
                }
            }
        }

        if(bestPartitionPoints == null)
            return null;

        List<List<ObservedRegion>> result = Lists.newArrayList();

        for(int j = 0; j < bestPartitionPoints.length - 1; j++)
        {
            result.add(Lists.newArrayList(members.subList(bestPartitionPoints[j], bestPartitionPoints[j + 1])));
        }

        return result;
    }

    static double deviationPenalty(final List<ObservedRegion> subsegmentMembers)
    {
        int m = subsegmentMembers.size();
        double mean = subsegmentMembers.stream().mapToDouble(x -> x.observedTumorRatio()).average().orElse(0);

        double total = 0;

        for(int i = 0; i < m; i++)
        {
            double ratio = subsegmentMembers.get(i).observedTumorRatio();
            double prevDiff = (i > 0) ? abs(ratio - subsegmentMembers.get(i - 1).observedTumorRatio()) : 0;
            double nextDiff = (i < m - 1) ? abs(ratio - subsegmentMembers.get(i + 1).observedTumorRatio()) : 0;
            double meanDiff = abs(ratio - mean);

            double deviation = max(prevDiff, max(nextDiff, meanDiff));
            total += deviation * deviation;
        }

        return total;
    }

    private static ObservedRegion aggregateRegions(final List<ObservedRegion> regions)
    {
        ObservedRegion first = regions.get(0);
        ObservedRegion last = regions.get(regions.size() - 1);

        boolean useBaf = regions.stream().anyMatch(x -> x.bafCount() >= RESEG_BAF_WEIGHT_THRESHOLD);

        int bafCountsTotal = 0;
        int depthWindowCountsTotal = 0;
        double weightedObsBafTotal = 0;
        double weightedObsTumorRatioTotal = 0;

        for(ObservedRegion region : regions)
        {
            int count = useBaf ? region.bafCount() : region.depthWindowCount();

            depthWindowCountsTotal += region.depthWindowCount();
            bafCountsTotal += region.bafCount();

            weightedObsBafTotal += region.observedBAF() * count;
            weightedObsTumorRatioTotal += region.observedTumorRatio() * count;
        }

        int countsTotal = useBaf ? bafCountsTotal : depthWindowCountsTotal;
        double calcObservedBAF = weightedObsBafTotal / countsTotal;
        double calObservedTumorRatio = weightedObsTumorRatioTotal / countsTotal;

        ObservedRegion newRegion = fromOther(first);
        newRegion.setEnd(last.end());
        newRegion.setBafCount(bafCountsTotal);
        newRegion.setObservedBAF(calcObservedBAF);
        newRegion.setDepthWindowCount(depthWindowCountsTotal);
        newRegion.setObservedTumorRatio(calObservedTumorRatio);

        return newRegion;
    }

    private static List<int[]> combinationsOf(int rangeSize, int k)
    {
        List<int[]> result = new ArrayList<>();
        combinationsHelper(result, new int[k], 0, 1, rangeSize, k);
        return result;
    }

    private static void combinationsHelper(final List<int[]> result, final int[] combo, int comboIndex, int start, int rangeSize, int k)
    {
        // finds all possible splits in segments
        if(comboIndex == k)
        {
            result.add(combo.clone());
            return;
        }

        for(int value = start; value <= rangeSize - (k - comboIndex - 1); value++)
        {
            combo[comboIndex] = value;
            combinationsHelper(result, combo, comboIndex + 1, value + 1, rangeSize, k);
        }
    }

    private static int[] pad(final int[] cuts, int n)
    {
        int[] padded = new int[cuts.length + 2];
        padded[0] = 0;
        arraycopy(cuts, 0, padded, 1, cuts.length);
        padded[padded.length - 1] = n;
        return padded;
    }
}
