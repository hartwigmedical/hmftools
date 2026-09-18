package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

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

    private static final String GROUP_A = "Group A";
    private static final String GROUP_H = "Group H";

    // Contigs are grouped by name prefix: "h" -> Group H, everything else -> Group A.
    private static final ViralReference REFERENCE = reference("v1", "v2", "v3", "h1", "h2");

    // Two contigs with similar votes and no challenge between them: the one with more votes leads, the other is its twin.
    @Test
    public void testResolvedTwins()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));

        List<OncologyGroupSelection> selections = select(stats, new Margins().build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, "v2"));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
        assertEquals(OncologyGroupResolution.RESOLVED, group(selections, GROUP_A).resolution());
        assertEquals(REFERENCE.contig("v1"), group(selections, GROUP_A).representative());
    }

    // The leader decisively challenges the other contig: the other becomes secondary and the leader is still chosen.
    @Test
    public void testResolvedWithSecondary()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));

        List<OncologyGroupSelection> selections = select(stats, new Margins().challenge("v1", "v2", 40).build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigRole.SECONDARY, role(selections, "v2"));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
        assertEquals(List.of(2), candidate(selections, "v1").challengesRanks());
        assertEquals(List.of(1), candidate(selections, "v2").challengedByRanks());
    }

    // Two contigs challenge each other, so neither can be chosen: unresolved, no representative.
    @Test
    public void testUnresolvedMutual()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 95));
        PairwiseMargins margins = new Margins().challenge("v1", "v2", 40).challenge("v2", "v1", 40).build();

        List<OncologyGroupSelection> selections = select(stats, margins);

        assertEquals(OncologyGroupOutcome.MUTUAL, group(selections, GROUP_A).outcome());
        assertEquals(OncologyGroupResolution.UNRESOLVED, group(selections, GROUP_A).resolution());
        assertNull(group(selections, GROUP_A).representative());
    }

    // Three contigs challenge in a loop (v1 beats v2 beats v3 beats v1), so none is unchallenged: unresolved as a cycle.
    @Test
    public void testUnresolvedCycle()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v2", 98), present("v3", 96));
        PairwiseMargins margins = new Margins()
                .challenge("v1", "v2", 40).challenge("v2", "v3", 40).challenge("v3", "v1", 40).build();

        List<OncologyGroupSelection> selections = select(stats, margins);

        assertEquals(OncologyGroupOutcome.CYCLE, group(selections, GROUP_A).outcome());
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig decisively challenges the leader: unresolved, and the leader is left contested.
    @Test
    public void testUnresolvedMinorChallenger()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v3", 10));

        List<OncologyGroupSelection> selections = select(stats, new Margins().challenge("v3", "v1", 40).build());

        assertEquals(OncologyGroupOutcome.MINOR_RIVAL, group(selections, GROUP_A).outcome());
        assertEquals(ContigRole.CONTESTED, role(selections, "v1"));
        assertEquals(ContigRole.MINOR_CHALLENGER, role(selections, "v3"));
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig that challenges nobody: the leader is still chosen and the low-vote contig is minor.
    @Test
    public void testResolvedWithMinorBystander()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100), present("v3", 10));

        List<OncologyGroupSelection> selections = select(stats, new Margins().build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigRole.MINOR, role(selections, "v3"));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
    }

    // A single candidate resolves trivially.
    @Test
    public void testSoleContig()
    {
        Map<String, ContigStats> stats = statsMap(present("v1", 100));

        List<OncologyGroupSelection> selections = select(stats, new Margins().build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(OncologyGroupOutcome.ONE_CANDIDATE, group(selections, GROUP_A).outcome());
    }

    // A contig below the relaxed coverage floor is dropped; a group where no contig meets the coverage minimum keeps nothing.
    @Test
    public void testPrefilteredContigsTakeNoPart()
    {
        Map<String, ContigStats> stats = statsMap(
                present("v1", 100),
                stats("v2", 0.05, 100),   // in a covered group, but itself below the relaxed floor
                stats("h1", 0.05, 100),   // no contig in this group meets the coverage minimum
                stats("h2", 0.04, 100));

        List<OncologyGroupSelection> selections = select(stats, new Margins().build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigFilterStatus.LOW_COVERAGE, contig(selections, "v2").filterStatus());
        assertNull(contig(selections, "v2").candidate());
        assertEquals(OncologyGroupOutcome.NO_CANDIDATES, group(selections, GROUP_H).outcome());
        assertEquals(OncologyGroupResolution.NO_CANDIDATES, group(selections, GROUP_H).resolution());
        assertNull(group(selections, GROUP_H).representative());
    }

    // A contig with good coverage but almost no votes is dropped by the vote-density floor.
    @Test
    public void testVoteDensityPrefilterDropsLowVoteContig()
    {
        Map<String, ContigStats> stats = statsMap(
                present("v1", 100),
                stats("v2", 0.5, 0.1));   // good coverage but almost no votes

        List<OncologyGroupSelection> selections = select(stats, new Margins().build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigFilterStatus.LOW_VOTE_DENSITY, contig(selections, "v2").filterStatus());
        assertNull(contig(selections, "v2").candidate());
    }

    private static List<OncologyGroupSelection> select(Map<String, ContigStats> stats, PairwiseMargins margins)
    {
        return new RepresentativeSelector().select(stats.values(), margins, groupReadCounts(stats), MEAN_READ_LENGTH);
    }

    // Each group is given as many reads as its contigs have votes, so a contig's votes read as its share of the group.
    private static Map<String, Integer> groupReadCounts(Map<String, ContigStats> stats)
    {
        Map<String, Integer> readCounts = new HashMap<>();
        stats.values().forEach(stat -> readCounts.merge(
                stat.contig().oncologyGroup(), (int) Math.round(stat.readVotes()), Integer::sum));
        return readCounts;
    }

    private static OncologyGroupSelection group(List<OncologyGroupSelection> selections, String oncologyGroup)
    {
        return selections.stream()
                .filter(selection -> selection.oncologyGroup().equals(oncologyGroup))
                .findFirst()
                .orElseThrow();
    }

    private static ContigSelectionResult contig(List<OncologyGroupSelection> selections, String contig)
    {
        return selections.stream()
                .flatMap(selection -> selection.contigs().stream())
                .filter(result -> result.contig().name().equals(contig))
                .findFirst()
                .orElseThrow();
    }

    private static CandidateSelectionResult candidate(List<OncologyGroupSelection> selections, String contig)
    {
        return contig(selections, contig).candidate();
    }

    private static ContigRole role(List<OncologyGroupSelection> selections, String contig)
    {
        return candidate(selections, contig).role();
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
            String group = name.startsWith("h") ? GROUP_H : GROUP_A;
            contigs.add(new ViralContig(name, LENGTH, "Virus " + name, group));
            records.add(new SAMSequenceRecord(name, LENGTH));
        }
        return new ViralReference(contigs, new SAMSequenceDictionary(records));
    }

    // Builds pairwise margins with controlled challenge edges: each challenge places its reads at a margin above the
    // selection margin, so they count as decisive.
    private static class Margins
    {
        private final Map<ContigPair, NavigableMap<Integer, Integer>> mMargins = new HashMap<>();
        private final Map<ContigPair, Integer> mShared = new HashMap<>();

        private Margins challenge(String subject, String opponent, int reads)
        {
            ContigPair pair = new ContigPair(REFERENCE.contig(subject), REFERENCE.contig(opponent));
            mMargins.computeIfAbsent(pair, key -> new TreeMap<>()).merge(10, reads, Integer::sum);
            mShared.merge(pair, reads, Integer::sum);
            mShared.merge(new ContigPair(REFERENCE.contig(opponent), REFERENCE.contig(subject)), reads, Integer::sum);
            return this;
        }

        private PairwiseMargins build()
        {
            return new PairwiseMargins(mMargins, mShared);
        }
    }
}
