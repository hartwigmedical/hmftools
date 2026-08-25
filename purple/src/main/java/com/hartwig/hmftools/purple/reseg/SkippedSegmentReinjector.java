package com.hartwig.hmftools.purple.reseg;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.purple.region.ObservedRegion.fromOther;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class SkippedSegmentReinjector
{
    public static List<ObservedRegion> reinjectSkippedRegions(final List<ObservedRegion> resegmentedRegions, final Supersegment supersegment)
    {
        if(supersegment.SkippableMembers.isEmpty())
        {
            if(resegmentedRegions.size() == supersegment.BothNoneMembers.size())
                return resegmentedRegions;

            if(resegmentedRegions.size() != 1)
                return List.of(formBlendedRegion(supersegment.BothNoneMembers, resegmentedRegions.get(0)));
        }

        List<ObservedRegion> result = Lists.newArrayList();

        for(ObservedRegion subsegment : resegmentedRegions)
        {
            List<ObservedRegion> skipsToInject = supersegment.SkippableMembers.stream()
                    .filter(skip -> skip.start() > subsegment.start() && skip.end() < subsegment.end())
                    .sorted(Comparator.comparingInt(s -> s.start()))
                    .collect(Collectors.toList());

            if(skipsToInject.isEmpty())
            {
                result.add(buildReinjectedSegment(subsegment, supersegment, subsegment.start(), subsegment.end()));
                continue;
            }

            List<int[]> ranges = Lists.newArrayList();
            ranges.add(new int[] { subsegment.start(), skipsToInject.get(0).start() - 1 });

            for(int i = 0; i < skipsToInject.size() - 1; i++)
            {
                ranges.add(new int[] { skipsToInject.get(i).end() + 1, skipsToInject.get(i + 1).start() - 1 });
            }

            ranges.add(new int[] { skipsToInject.get(skipsToInject.size() - 1).end() + 1, subsegment.end() });

            for(int[] range : ranges)
            {
                result.add(buildReinjectedSegment(subsegment, supersegment, range[0], range[1]));
            }
        }

        return result;
    }

    private static ObservedRegion buildReinjectedSegment(
            final ObservedRegion subsegment, final Supersegment supersegment, int rangeStart, int rangeEnd)
    {
        List<ObservedRegion> relevantBothNone = supersegment.BothNoneMembers.stream()
                .filter(s -> s.start() >= rangeStart && s.end() <= rangeEnd)
                .collect(Collectors.toList());

        if(relevantBothNone.isEmpty())
        {
            throw new IllegalStateException(format("no bothNone members found for reinjected range(%d-%d)", rangeStart, rangeEnd));
        }

        if(relevantBothNone.size() == 1)
        {
            ObservedRegion newRegion = fromOther(relevantBothNone.get(0));
            newRegion.setObservedTumorRatio(subsegment.observedTumorRatio());
            newRegion.setObservedBAF(subsegment.observedBAF());
            return newRegion;
        }

        return formBlendedRegion(relevantBothNone, subsegment);
    }

    public static ObservedRegion formBlendedRegion(final List<ObservedRegion> regions, final ObservedRegion subsegment)
    {
        int bafTotal = 0;
        int depthWindowTotal = 0;
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

        newRegion.setObservedTumorRatio(subsegment.observedTumorRatio());
        newRegion.setObservedBAF(subsegment.observedBAF());

        newRegion.setBafCount(bafTotal);
        newRegion.setDepthWindowCount(depthWindowTotal);

        newRegion.setObservedNormalRatio(weightedObsNormRatioTotal / depthWindowTotal);
        newRegion.setUnnormalisedObservedNormalRatio(weightedUnnormalisedObsNormRatioTotal / depthWindowTotal);
        newRegion.setRatioSupport(hasRatioSupport);
        newRegion.setGcContent(gcRatioTotal / totalLength);

        return newRegion;
    }
}
