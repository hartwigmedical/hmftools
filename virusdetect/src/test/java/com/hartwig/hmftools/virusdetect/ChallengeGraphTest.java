package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import com.hartwig.hmftools.virusdetect.PairwiseMargins.ContigPair;

import org.junit.Test;

public class ChallengeGraphTest
{
    private static final ViralContig V1 = new ViralContig("v1", 1000, "Virus v1", new OncologyGroup("Group A"));
    private static final ViralContig V2 = new ViralContig("v2", 1000, "Virus v2", new OncologyGroup("Group A"));
    private static final ViralContig V3 = new ViralContig("v3", 1000, "Virus v3", new OncologyGroup("Group A"));

    // A tenth of the group's reads win by a decisive margin, which meets the challenge threshold of a tenth.
    @Test
    public void testChallengeAtThreshold()
    {
        List<ContigSupport> group = List.of(stats(V1, 60), stats(V2, 40));
        ChallengeGraph graph = ChallengeGraph.build(group, 100, margins(V1, V2, MIN_CHALLENGE_MARGIN, 10));

        assertTrue(graph.challenges(V1, V2));
        assertFalse(graph.challenges(V2, V1));
    }

    @Test
    public void testNoChallengeBelowThreshold()
    {
        List<ContigSupport> group = List.of(stats(V1, 60), stats(V2, 40));
        ChallengeGraph graph = ChallengeGraph.build(group, 100, margins(V1, V2, MIN_CHALLENGE_MARGIN, 9));

        assertFalse(graph.challenges(V1, V2));
    }

    // Plenty of winning reads, but none of them wins by enough bases to count.
    @Test
    public void testNoChallengeBelowMargin()
    {
        List<ContigSupport> group = List.of(stats(V1, 60), stats(V2, 40));
        ChallengeGraph graph = ChallengeGraph.build(group, 100, margins(V1, V2, MIN_CHALLENGE_MARGIN - 1, 50));

        assertFalse(graph.challenges(V1, V2));
    }

    // The fraction is over every read aligning to the group, including reads which only the contigs dropped before
    // selection carry. These 8 reads would be 1/10 of the two candidates' own reads, but are a 1/12 of the group's.
    @Test
    public void testFractionSpansEveryGroupRead()
    {
        List<ContigSupport> candidates = List.of(stats(V1, 40), stats(V2, 40));
        ChallengeGraph graph = ChallengeGraph.build(candidates, 96, margins(V1, V2, MIN_CHALLENGE_MARGIN, 8));

        assertFalse(graph.challenges(V1, V2));
    }

    // Abundant contigs are those whose votes are near the highest.
    @Test
    public void testComparableContigs()
    {
        List<ContigSupport> group = List.of(stats(V1, 100), stats(V2, 95), stats(V3, 50));
        ChallengeGraph graph = ChallengeGraph.build(group, 245, margins(V1, V2, MIN_CHALLENGE_MARGIN, 0));

        assertEquals(Set.of(V1, V2), graph.comparable());
    }

    private static PairwiseMargins margins(ViralContig subject, ViralContig opponent, int margin, int reads)
    {
        ContigPair pair = new ContigPair(subject, opponent);
        Map<ContigPair, NavigableMap<Integer, Integer>> marginCounts = new HashMap<>();
        marginCounts.put(pair, new TreeMap<>(Map.of(margin, reads)));
        Map<ContigPair, Integer> sharedReads = new HashMap<>();
        sharedReads.put(pair, reads);
        return new PairwiseMargins(marginCounts, sharedReads);
    }

    private static ContigSupport stats(ViralContig contig, double votes)
    {
        SummaryStats dummy = SummaryStats.from(new int[] { 1 });
        return new ContigSupport(contig, 100, 0, dummy, 0, 500, dummy, dummy, votes);
    }
}
