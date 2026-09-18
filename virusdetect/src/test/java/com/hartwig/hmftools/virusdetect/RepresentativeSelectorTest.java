package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import com.hartwig.hmftools.virusdetect.PairwiseMargins.ContigPair;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class RepresentativeSelectorTest
{
    // All contigs are length 1000 with mean read length 150, so the vote-density floor is
    // 0.10 * 1000 / 150 ~= 0.67 votes, which every contig here clears unless it has almost no votes.
    private static final int LENGTH = 1000;
    private static final double MEAN_READ_LENGTH = 150.0;

    // Contigs are grouped by name prefix: "h" -> Group H, everything else -> Group A.
    private static final ViralReference REFERENCE = reference("v1", "v2", "v3", "h1", "h2");

    // Two contigs with similar votes and no challenge between them: the one with more votes leads, the other is its twin.
    @Test
    public void testResolvedTwins()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));
        PairwiseMargins pairwise = new Pairwise().build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(ContigRole.REPRESENTATIVE_TWIN, byContig.get("v2").role());
        assertEquals(OncologyGroupSubOutcome.RESOLVED_CANDIDATES, byContig.get("v1").oncologyGroupSubOutcome());
        assertEquals(OncologyGroupOutcome.RESOLVED, byContig.get("v1").oncologyGroupOutcome());
    }

    // The leader decisively challenges the other contig: the other becomes secondary and the leader is still chosen.
    @Test
    public void testResolvedWithSecondary()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));
        PairwiseMargins pairwise = new Pairwise().challenge("v1", "v2", 40).build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(ContigRole.SECONDARY, byContig.get("v2").role());
        assertEquals(OncologyGroupSubOutcome.RESOLVED_CANDIDATES, byContig.get("v1").oncologyGroupSubOutcome());
        assertEquals(List.of(2), byContig.get("v1").challengesRanks());
        assertEquals(List.of(1), byContig.get("v2").challengedByRanks());
    }

    // Two contigs challenge each other, so neither can be chosen: unresolved, no representative.
    @Test
    public void testUnresolvedMutual()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));
        PairwiseMargins pairwise = new Pairwise().challenge("v1", "v2", 40).challenge("v2", "v1", 40).build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(OncologyGroupSubOutcome.MUTUAL, byContig.get("v1").oncologyGroupSubOutcome());
        assertEquals(OncologyGroupOutcome.UNRESOLVED, byContig.get("v1").oncologyGroupOutcome());
        assertNoRepresentative(byContig.values());
    }

    // Three contigs challenge in a loop (v1 beats v2 beats v3 beats v1), so none is unchallenged: unresolved as a cycle.
    @Test
    public void testUnresolvedCycle()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 98), present("v3", 96));
        PairwiseMargins pairwise = new Pairwise()
                .challenge("v1", "v2", 40).challenge("v2", "v3", 40).challenge("v3", "v1", 40).build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(OncologyGroupSubOutcome.CYCLE, byContig.get("v1").oncologyGroupSubOutcome());
        assertNoRepresentative(byContig.values());
    }

    // A low-vote contig decisively challenges the leader: unresolved, and the leader is left contested.
    @Test
    public void testUnresolvedMinorChallenger()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v3", 10));
        PairwiseMargins pairwise = new Pairwise().challenge("v3", "v1", 40).build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(OncologyGroupSubOutcome.MINOR_RIVAL, byContig.get("v1").oncologyGroupSubOutcome());
        assertEquals(ContigRole.CONTESTED, byContig.get("v1").role());
        assertEquals(ContigRole.MINOR_CHALLENGER, byContig.get("v3").role());
        assertNoRepresentative(byContig.values());
    }

    // A low-vote contig that challenges nobody: the leader is still chosen and the low-vote contig is minor.
    @Test
    public void testResolvedWithMinorBystander()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v3", 10));
        PairwiseMargins pairwise = new Pairwise().build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(ContigRole.MINOR, byContig.get("v3").role());
        assertEquals(OncologyGroupSubOutcome.RESOLVED_CANDIDATES, byContig.get("v1").oncologyGroupSubOutcome());
    }

    // A single candidate resolves trivially.
    @Test
    public void testSoleContig()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100));
        PairwiseMargins pairwise = new Pairwise().build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(OncologyGroupSubOutcome.ONE_CANDIDATE, byContig.get("v1").oncologyGroupSubOutcome());
    }

    // A contig below the relaxed coverage floor is dropped; a group where no contig meets the coverage minimum keeps nothing.
    @Test
    public void testCoveragePrefilter()
    {
        Map<String, ContigStats> stats = statsMap(
                present("v1", 100),
                stats("v2", 0.05, 100),   // in a covered group, but itself below the relaxed floor
                stats("h1", 0.05, 100),   // no contig in this group meets the coverage minimum
                stats("h2", 0.04, 100));
        PairwiseMargins pairwise = new Pairwise().build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(ContigFilterStatus.LOW_COVERAGE, byContig.get("v2").filterStatus());
        assertEquals(ContigFilterStatus.LOW_COVERAGE, byContig.get("h1").filterStatus());
        assertEquals(ContigFilterStatus.LOW_COVERAGE, byContig.get("h2").filterStatus());
    }

    // A contig with good coverage but almost no votes is dropped by the vote-density floor.
    @Test
    public void testVoteDensityPrefilterDropsLowVoteContig()
    {
        Map<String, ContigStats> stats = statsMap(
                present("v1", 100),
                stats("v2", 0.5, 0.1));   // good coverage but almost no votes
        PairwiseMargins pairwise = new Pairwise().build();

        Map<String, ContigClassification> byContig = classify(stats, pairwise);

        assertEquals(ContigRole.REPRESENTATIVE, byContig.get("v1").role());
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, byContig.get("v2").filterStatus());
        // Still ranked among the covered contigs, but gets no role
        assertEquals(Integer.valueOf(2), byContig.get("v2").votesRank());
    }

    private static void assertNoRepresentative(Iterable<ContigClassification> classifications)
    {
        for(ContigClassification classification : classifications)
        {
            assertTrue("expected no representative", classification.role() != ContigRole.REPRESENTATIVE);
        }
    }

    private Map<String, ContigClassification> classify(Map<String, ContigStats> stats, PairwiseMargins pairwise)
    {
        RepresentativeSelectionResult result = new RepresentativeSelector().classify(stats.values(), pairwise);
        Map<String, ContigClassification> byContig = new HashMap<>();
        result.classifications().forEach(classification -> byContig.put(classification.contig().name(), classification));
        return byContig;
    }

    private static Map<String, ContigStats> statsMap(ContigStats... stats)
    {
        Map<String, ContigStats> map = new HashMap<>();
        for(ContigStats stat : stats)
        {
            map.put(stat.contig().name(), stat);
        }
        return map;
    }

    // A contig in Group A with good coverage and the given read votes.
    private static ContigStats present(String contig, double votes)
    {
        return stats(contig, 0.5, votes);
    }

    private static ContigStats stats(String contig, double coverage, double votes)
    {
        int coveredBases = (int) Math.round(coverage * LENGTH);
        SummaryStats dummy = SummaryStats.from(new int[] { 1 });
        return new ContigStats(REFERENCE.contig(contig), 100, 0, dummy, 0, coveredBases, dummy, dummy, votes);
    }

    private static ViralReference reference(String... contigNames)
    {
        List<ViralContig> contigs = new ArrayList<>();
        List<SAMSequenceRecord> records = new ArrayList<>();
        for(String name : contigNames)
        {
            String group = name.startsWith("h") ? "Group H" : "Group A";
            contigs.add(new ViralContig(name, LENGTH, "Virus " + name, group));
            records.add(new SAMSequenceRecord(name, LENGTH));
        }
        return new ViralReference(contigs, new SAMSequenceDictionary(records));
    }

    // Builds a PairwiseMargins with controlled challenge edges: each challenge places its reads at a margin above the
    // selection margin, so they count as decisive.
    private static class Pairwise
    {
        private final Map<ContigPair, NavigableMap<Integer, Integer>> mMargins = new HashMap<>();
        private final Map<ContigPair, Integer> mShared = new HashMap<>();

        private Pairwise challenge(String subject, String opponent, int reads)
        {
            ContigPair pair = new ContigPair(REFERENCE.contig(subject), REFERENCE.contig(opponent));
            mMargins.computeIfAbsent(pair, key -> new TreeMap<>()).merge(10, reads, Integer::sum);
            mShared.merge(pair, reads, Integer::sum);
            mShared.merge(new ContigPair(REFERENCE.contig(opponent), REFERENCE.contig(subject)), reads, Integer::sum);
            return this;
        }

        private PairwiseMargins build()
        {
            return new PairwiseMargins(mMargins, mShared, MEAN_READ_LENGTH);
        }
    }
}
