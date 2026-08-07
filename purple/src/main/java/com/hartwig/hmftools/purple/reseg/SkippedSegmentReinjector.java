package com.hartwig.hmftools.purple.reseg;

import static java.lang.String.format;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public final class SkippedSegmentReinjector
{
    private SkippedSegmentReinjector() {}

    public static List<ObservedRegion> reinject(final List<ObservedRegion> step5Subsegments, final Supersegment supersegment)
    {
        if(supersegment.SkippableMembers.isEmpty())
            return step5Subsegments;

        List<ObservedRegion> result = new ArrayList<>();

        for(ObservedRegion subsegment : step5Subsegments)
        {
            List<ObservedRegion> skipsToInject = supersegment.SkippableMembers.stream()
                    .filter(skip -> skip.start() > subsegment.start() && skip.end() < subsegment.end())
                    .sorted(Comparator.comparingInt(s -> s.start()))
                    .collect(Collectors.toList());

            if(skipsToInject.isEmpty())
            {
                result.add(subsegment);
                continue;
            }

            List<int[]> ranges = new ArrayList<>();
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
            throw new IllegalStateException(format(
                    "no bothNone members found for reinjected range(%d-%d)", rangeStart, rangeEnd));
        }

        int minStart = relevantBothNone.stream().mapToInt(x -> x.start()).min().getAsInt();
        int maxEnd = relevantBothNone.stream().mapToInt(x -> x.end()).max().getAsInt();

        if(minStart != rangeStart || maxEnd != rangeEnd)
        {
            throw new IllegalStateException(format(
                    "reinjected range(%d-%d) bounds mismatch: bothNone members span(%d-%d)",
                    rangeStart, rangeEnd, minStart, maxEnd));
        }

        int bafCount = relevantBothNone.stream().mapToInt(x -> x.bafCount()).sum();
        int depthWindowCount = relevantBothNone.stream().mapToInt(x -> x.depthWindowCount()).sum();

        ObservedRegion newRegion = ObservedRegion.fromOther(subsegment);
        newRegion.setStart(rangeStart);
        newRegion.setEnd(rangeEnd);
        newRegion.setBafCount(bafCount);
        newRegion.setDepthWindowCount(depthWindowCount);

        return newRegion;
    }
}
