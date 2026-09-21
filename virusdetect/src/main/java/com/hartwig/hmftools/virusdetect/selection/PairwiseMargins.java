package com.hartwig.hmftools.virusdetect.selection;

import static java.util.stream.Collectors.groupingBy;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import com.hartwig.hmftools.virusdetect.detection.read_align.ReadAlignments;
import com.hartwig.hmftools.virusdetect.detection.read_align.ViralReadAlignments;
import com.hartwig.hmftools.virusdetect.reference.OncologyGroup;
import com.hartwig.hmftools.virusdetect.reference.ViralContig;

// Per ordered within-group contig pair, how decisively the subject fits shared reads better than the opponent.
// Quantifies the distribution of positive winning alignment margins and the read count.
// This is a key input to the representative contig selection scheme, used to determine if there is significant support
// for a rival virus strain.
public class PairwiseMargins
{
    // Subject -> opponent -> (winning margin in bases -> read count)
    private final Map<ContigPair, NavigableMap<Integer, Integer>> mMarginCounts;
    private final Map<ContigPair, Integer> mSharedReads;

    public PairwiseMargins(Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts, Map<ContigPair, Integer> sharedReads)
    {
        mMarginCounts = marginCounts;
        mSharedReads = sharedReads;
    }

    public static PairwiseMargins from(ViralReadAlignments viralAlignments)
    {
        Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts = new HashMap<>();
        Map<ContigPair, Integer> sharedReads = new HashMap<>();

        viralAlignments.reads().forEach(read -> accumulateRead(read, marginCounts, sharedReads));

        return new PairwiseMargins(marginCounts, sharedReads);
    }

    // Reads whose best alignment fits the subject at least minMargin divergent bases better than the opponent.
    public int readsWinningBy(ViralContig subject, ViralContig opponent, int minMargin)
    {
        NavigableMap<Integer, Integer> margins = mMarginCounts.get(new ContigPair(subject, opponent));
        if(margins == null)
        {
            return 0;
        }
        else
        {
            return margins.tailMap(minMargin, true).values().stream().mapToInt(Integer::intValue).sum();
        }
    }

    // Reads aligned to both contigs (regardless of which one they favour).
    public int sharedReads(ViralContig subject, ViralContig opponent)
    {
        return mSharedReads.getOrDefault(new ContigPair(subject, opponent), 0);
    }

    public Set<ContigPair> pairs()
    {
        return mSharedReads.keySet();
    }

    // Considers one read at a time.
    // Pairs up the contigs it aligns to within each oncology group.
    // For each ordered pair, records that they share the read, and the subject's winning margin over the opponent (if any).
    private static void accumulateRead(
            ReadAlignments read,
            Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts, Map<ContigPair, Integer> sharedReads)
    {
        Map<OncologyGroup, List<ViralContig>> contigsByOncologyGroup = read.hits().keySet().stream()
                .collect(groupingBy(ViralContig::oncologyGroup));

        for(List<ViralContig> oncologyGroupContigs : contigsByOncologyGroup.values())
        {
            for(ViralContig subject : oncologyGroupContigs)
            {
                for(ViralContig opponent : oncologyGroupContigs)
                {
                    if(!subject.equals(opponent))
                    {
                        ContigPair pair = new ContigPair(subject, opponent);
                        sharedReads.merge(pair, 1, Integer::sum);

                        int margin = read.hits().get(opponent).divergence() - read.hits().get(subject).divergence();
                        if(margin > 0)
                        {
                            marginCounts.computeIfAbsent(pair, key -> new TreeMap<>()).merge(margin, 1, Integer::sum);
                        }
                    }
                }
            }
        }
    }

    public record ContigPair(
            ViralContig subject,
            ViralContig opponent
    )
    {
    }
}
