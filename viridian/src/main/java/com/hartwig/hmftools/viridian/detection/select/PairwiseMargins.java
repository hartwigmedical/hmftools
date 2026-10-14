package com.hartwig.hmftools.viridian.detection.select;

import static java.util.stream.Collectors.groupingBy;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.common.ReadId;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Per ordered within-group virus genome pair, how decisively the subject fits shared reads better than the opponent.
// Quantifies the distribution of positive winning alignment margins and the read count.
// This is a key input to the representative selection scheme, used to determine if there is significant support
// for a rival virus strain.
public class PairwiseMargins
{
    // Subject -> opponent -> (winning margin in bases -> read count)
    private final Map<GenomePair, NavigableMap<Integer, Integer>> mMarginCounts;
    private final Map<GenomePair, Integer> mSharedReads;

    public PairwiseMargins(Map<GenomePair, NavigableMap<Integer, Integer>> marginCounts, Map<GenomePair, Integer> sharedReads)
    {
        mMarginCounts = marginCounts;
        mSharedReads = sharedReads;
    }

    public static PairwiseMargins from(Map<ReadId, Map<VirusGenome, ViralReadAlignment>> alignmentsByRead)
    {
        Map<GenomePair, NavigableMap<Integer, Integer>> marginCounts = new HashMap<>();
        Map<GenomePair, Integer> sharedReads = new HashMap<>();

        alignmentsByRead.values().forEach(read -> accumulateRead(read, marginCounts, sharedReads));

        return new PairwiseMargins(marginCounts, sharedReads);
    }

    // Reads whose best alignment fits the subject at least minMargin divergent bases better than the opponent.
    public int readsWinningBy(VirusGenome subject, VirusGenome opponent, int minMargin)
    {
        NavigableMap<Integer, Integer> margins = mMarginCounts.get(new GenomePair(subject, opponent));
        if(margins == null)
        {
            return 0;
        }
        else
        {
            return margins.tailMap(minMargin, true).values().stream().mapToInt(Integer::intValue).sum();
        }
    }

    // Reads aligned to both virus genomes (regardless of which one they favour).
    public int sharedReads(VirusGenome subject, VirusGenome opponent)
    {
        return mSharedReads.getOrDefault(new GenomePair(subject, opponent), 0);
    }

    public Set<GenomePair> pairs()
    {
        return mSharedReads.keySet();
    }

    // Considers one read at a time.
    // Pairs up the virus genomes it aligns to within each oncology group.
    // For each ordered pair, records that they share the read, and the subject's winning margin over the opponent (if any).
    private static void accumulateRead(
            Map<VirusGenome, ViralReadAlignment> readAlignments,
            Map<GenomePair, NavigableMap<Integer, Integer>> marginCounts, Map<GenomePair, Integer> sharedReads)
    {
        Map<OncologyGroup, List<VirusGenome>> genomesByOncologyGroup = readAlignments.keySet().stream()
                .collect(groupingBy(VirusGenome::oncologyGroup));

        for(List<VirusGenome> oncologyGroupGenomes : genomesByOncologyGroup.values())
        {
            for(VirusGenome subject : oncologyGroupGenomes)
            {
                for(VirusGenome opponent : oncologyGroupGenomes)
                {
                    if(!subject.equals(opponent))
                    {
                        GenomePair pair = new GenomePair(subject, opponent);
                        sharedReads.merge(pair, 1, Integer::sum);

                        int margin = readAlignments.get(opponent).divergence() - readAlignments.get(subject).divergence();
                        if(margin > 0)
                        {
                            marginCounts.computeIfAbsent(pair, k -> new TreeMap<>()).merge(margin, 1, Integer::sum);
                        }
                    }
                }
            }
        }
    }

    public record GenomePair(
            VirusGenome subject,
            VirusGenome opponent
    )
    {
    }
}
