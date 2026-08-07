package com.hartwig.hmftools.purple.reseg;

import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_BAF_WEIGHT_THRESHOLD;
import static com.hartwig.hmftools.purple.PurpleConstants.RESEG_MAX_SUPERSEGMENT_SIZE_FOR_BRUTE_FORCE;
import static com.hartwig.hmftools.purple.region.ObservedRegion.fromOther;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleFunction;
import com.hartwig.hmftools.purple.region.ObservedRegion;

/**
 * Step 5: resegments a supersegment's members into subsegments, trading off per-subsegment
 * observedTumorRatio deviation against the per-sample segmentation penalty.
 */
public final class SupersegmentResegmenter
{
    private SupersegmentResegmenter() {}

    public static List<ObservedRegion> resegment(final List<ObservedRegion> members, double segmentationPenalty)
    {
        if(members.size() > RESEG_MAX_SUPERSEGMENT_SIZE_FOR_BRUTE_FORCE)
            return splitLargeSupersegment(members, segmentationPenalty);

        List<List<ObservedRegion>> bestPartition = findBestPartition(members, segmentationPenalty);

        if(bestPartition == null)
            return List.of(aggregate(members));

        List<ObservedRegion> result = new ArrayList<>();

        for(List<ObservedRegion> subsegment : bestPartition)
        {
            result.add(aggregate(subsegment));
        }

        return result;
    }

    private static List<ObservedRegion> splitLargeSupersegment(final List<ObservedRegion> members, double segmentationPenalty)
    {
        double maxRatio = members.stream().mapToDouble(m -> m.observedTumorRatio()).max().orElse(0);
        double minRatio = members.stream().mapToDouble(m -> m.observedTumorRatio()).min().orElse(0);

        if(maxRatio - minRatio > Math.sqrt(segmentationPenalty / 2))
        {
            int cutoffIndex = findLargestAdjacentDiffIndex(members);

            List<ObservedRegion> left = new ArrayList<>(members.subList(0, cutoffIndex));
            List<ObservedRegion> right = new ArrayList<>(members.subList(cutoffIndex, members.size()));

            List<ObservedRegion> result = new ArrayList<>();
            result.addAll(resegment(left, segmentationPenalty));
            result.addAll(resegment(right, segmentationPenalty));
            return result;
        }

        return List.of(aggregate(members));
    }

    private static int findLargestAdjacentDiffIndex(final List<ObservedRegion> members)
    {
        int bestIndex = 1;
        double bestDiff = -1;

        for(int i = 1; i < members.size(); i++)
        {
            double diff = Math.abs(members.get(i).observedTumorRatio() - members.get(i - 1).observedTumorRatio());

            if(diff >= bestDiff)
            {
                bestDiff = diff;
                bestIndex = i;
            }
        }

        return bestIndex;
    }

    // returns null if no partition beats the unpartitioned (0-cut) baseline - caller should aggregate members as one
    private static List<List<ObservedRegion>> findBestPartition(final List<ObservedRegion> members, double segmentationPenalty)
    {
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

        List<List<ObservedRegion>> result = new ArrayList<>();

        for(int j = 0; j < bestPartitionPoints.length - 1; j++)
        {
            result.add(new ArrayList<>(members.subList(bestPartitionPoints[j], bestPartitionPoints[j + 1])));
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
            double prevDiff = (i > 0) ? Math.abs(ratio - subsegmentMembers.get(i - 1).observedTumorRatio()) : 0;
            double nextDiff = (i < m - 1) ? Math.abs(ratio - subsegmentMembers.get(i + 1).observedTumorRatio()) : 0;
            double meanDiff = Math.abs(ratio - mean);

            double deviation = Math.max(prevDiff, Math.max(nextDiff, meanDiff));
            total += deviation * deviation;
        }

        return total;
    }

    static ObservedRegion aggregate(final List<ObservedRegion> members)
    {
        ObservedRegion first = members.get(0);
        ObservedRegion last = members.get(members.size() - 1);

        boolean useBaf = members.stream().mapToInt(m -> m.bafCount()).max().orElse(0) >= RESEG_BAF_WEIGHT_THRESHOLD;

        double weightSum = useBaf
                ? members.stream().mapToDouble(m -> m.bafCount()).sum()
                : members.stream().mapToDouble(m -> m.depthWindowCount()).sum();

        double observedBAF;
        double observedTumorRatio;

        if(weightSum > 0)
        {
            observedBAF = weightedSum(members, useBaf, m -> m.observedBAF()) / weightSum;
            observedTumorRatio = weightedSum(members, useBaf, m -> m.observedTumorRatio()) / weightSum;
        }
        else
        {
            observedBAF = members.stream().mapToDouble(m -> m.observedBAF()).average().orElse(0);
            observedTumorRatio = members.stream().mapToDouble(m -> m.observedTumorRatio()).average().orElse(0);
        }

        int bafCount = members.stream().mapToInt(m -> m.bafCount()).sum();
        int depthWindowCount = members.stream().mapToInt(m -> m.depthWindowCount()).sum();

        ObservedRegion newRegion = fromOther(first);
        newRegion.setEnd(last.end());
        newRegion.setBafCount(bafCount);
        newRegion.setObservedBAF(observedBAF);
        newRegion.setDepthWindowCount(depthWindowCount);
        newRegion.setObservedTumorRatio(observedTumorRatio);

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

    // all k-combinations of {1, ..., rangeSize}, in ascending lexicographic order (matches itertools.combinations)
    private static List<int[]> combinationsOf(int rangeSize, int k)
    {
        List<int[]> result = new ArrayList<>();
        combinationsHelper(result, new int[k], 0, 1, rangeSize, k);
        return result;
    }

    private static void combinationsHelper(final List<int[]> result, final int[] combo, int comboIndex, int start, int rangeSize, int k)
    {
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
        System.arraycopy(cuts, 0, padded, 1, cuts.length);
        padded[padded.length - 1] = n;
        return padded;
    }
}
