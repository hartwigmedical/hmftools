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
            return List.of(formNewRegion(members));

        List<ObservedRegion> result = new ArrayList<>();

        for(List<ObservedRegion> subsegment : bestPartition)
        {
            result.add(formNewRegion(subsegment));
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

        return List.of(formNewRegion(members));
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

    private static ObservedRegion formNewRegion(final List<ObservedRegion> regions)
    {
        int bafTotal = 0;
        int depthWindowTotal = 0;
        double weightedObsBAFTotal = 0;
        double weightedObsNormRatioTotal = 0;
        double weightedUnnormalisedObsNormRatioTotal = 0;
        double gcRatioTotal = 0;
        int totalLength = 0;

        boolean hasRatioSupport = false;
        int minStart = -1;
        int maxStart = 0;

        for(ObservedRegion region : regions)
        {
            bafTotal += region.bafCount();
            weightedObsBAFTotal += region.bafCount() * region.observedBAF();

            depthWindowTotal += region.depthWindowCount();
            weightedObsNormRatioTotal += region.depthWindowCount() * region.observedNormalRatio();

            weightedUnnormalisedObsNormRatioTotal += region.depthWindowCount() * region.unnormalisedObservedNormalRatio();

            hasRatioSupport |= region.ratioSupport();

            int regionLength = region.end() - region.start() + 1;
            totalLength += regionLength;
            gcRatioTotal += region.gcContent() * regionLength;

            minStart = minStart < 0 ? region.minStart() : min(minStart, region.minStart());
            maxStart = max(maxStart, region.maxStart());
        }

        ObservedRegion first = regions.get(0);
        ObservedRegion last = regions.get(regions.size() - 1);

        ObservedRegion newRegion = fromOther(first);
        newRegion.setEnd(last.end());
        newRegion.setMinStart(minStart);
        newRegion.setMaxStart(maxStart);

        newRegion.setBafCount(bafTotal);
        newRegion.setDepthWindowCount(depthWindowTotal);

        newRegion.setObservedBAF(weightedObsBAFTotal / bafTotal);
        newRegion.setObservedNormalRatio(weightedObsNormRatioTotal / depthWindowTotal);
        newRegion.setUnnormalisedObservedNormalRatio(weightedUnnormalisedObsNormRatioTotal / depthWindowTotal);
        newRegion.setRatioSupport(hasRatioSupport);
        newRegion.setGcContent(gcRatioTotal / totalLength);

        return newRegion;
    }

    private static double weightedSum(final List<ObservedRegion> members, boolean useBaf, final ToDoubleFunction<ObservedRegion> fieldFn)
    {
        double sum = 0;

        for(ObservedRegion m : members)
        {
            double weight = useBaf ? m.bafCount() : m.depthWindowCount();
            sum += weight * fieldFn.applyAsDouble(m);
        }

        return sum;
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
