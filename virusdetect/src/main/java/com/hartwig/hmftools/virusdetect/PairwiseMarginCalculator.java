package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.groupingBy;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import com.hartwig.hmftools.virusdetect.PairwiseMargins.ContigPair;

// For each pair (A, B) of viral contigs in an oncology group, looks at reads which align to both A and B.
// Calculates the "margin" of edit distance between the alignment to A and B.
// This is used to detect when the sample's virus strain is not well represented by any single contig in our reference.
// so a subset of reads will preferentially align to one contig with large margin.
public class PairwiseMarginCalculator
{
    public PairwiseMargins compute(ViralAlignments viralAlignments, ViralReference reference)
    {
        List<ViralAlignment> withinContig = viralAlignments.alignments().stream()
                .filter(alignment -> !alignment.clipsOverContigEnd(reference.contig(alignment.contig()).length()))
                .toList();

        Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts = new HashMap<>();
        Map<ContigPair, Integer> sharedReads = new HashMap<>();

        withinContig.stream().collect(groupingBy(ViralAlignment::readName)).values()
                .forEach(readAlignments -> accumulateRead(readAlignments, reference, marginCounts, sharedReads));

        return new PairwiseMargins(marginCounts, sharedReads, viralAlignments.meanReadLength());
    }

    // Considers one read at a time.
    // Pairs up the contigs it aligns to within each oncology group.
    // For each ordered pair, records that they share the read, and the subject's winning margin over the opponent (if any).
    private static void accumulateRead(
            List<ViralAlignment> readAlignments, ViralReference reference,
            Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts, Map<ContigPair, Integer> sharedReads)
    {
        Map<String, Integer> bestDivergenceByContig = new HashMap<>();
        readAlignments.forEach(alignment -> bestDivergenceByContig.merge(alignment.contig(), alignment.divergence(), Math::min));

        Map<String, List<String>> contigsByOncologyGroup = bestDivergenceByContig.keySet().stream()
                .collect(groupingBy(contig -> reference.contig(contig).oncologyGroup()));

        for(List<String> oncologyGroupContigs : contigsByOncologyGroup.values())
        {
            for(String subject : oncologyGroupContigs)
            {
                for(String opponent : oncologyGroupContigs)
                {
                    if(subject.equals(opponent))
                    {
                        continue;
                    }
                    ContigPair pair = new ContigPair(subject, opponent);
                    sharedReads.merge(pair, 1, Integer::sum);

                    int margin = bestDivergenceByContig.get(opponent) - bestDivergenceByContig.get(subject);
                    if(margin > 0)
                    {
                        marginCounts.computeIfAbsent(pair, key -> new TreeMap<>()).merge(margin, 1, Integer::sum);
                    }
                }
            }
        }
    }
}
