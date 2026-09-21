package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import com.hartwig.hmftools.virusdetect.PairwiseMargins.ContigPair;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class RepresentativeSelectorTest
{
    private static final int LENGTH = 1000;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    // Contigs are grouped by name prefix: "h" -> Group H, everything else -> Group A.
    private static final ViralReference REFERENCE = reference("v1", "v2", "v3", "h1", "h2");

    // Two contigs with similar votes and no challenge between them: the one with more votes leads, the other is its twin.
    @Test
    public void testResolvedTwins()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate("v1", 100), candidate("v2", 95)), noChallenges());

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
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate("v1", 100), candidate("v2", 95)), new Margins().challenge("v1", "v2", 40).build());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigRole.SECONDARY, role(selections, "v2"));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
        assertEquals(Set.of(REFERENCE.contig("v2")), result(selections, "v1").challenges());
        assertEquals(Set.of(REFERENCE.contig("v1")), result(selections, "v2").challengedBy());
    }

    // Two contigs challenge each other, so neither can be chosen: unresolved, no representative.
    @Test
    public void testUnresolvedMutual()
    {
        PairwiseMargins margins = new Margins().challenge("v1", "v2", 40).challenge("v2", "v1", 40).build();

        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate("v1", 100), candidate("v2", 95)), margins);

        assertEquals(OncologyGroupOutcome.MUTUAL, group(selections, GROUP_A).outcome());
        assertEquals(OncologyGroupResolution.UNRESOLVED, group(selections, GROUP_A).resolution());
        assertNull(group(selections, GROUP_A).representative());
    }

    // Three contigs challenge in a loop (v1 beats v2 beats v3 beats v1), so none is unchallenged: unresolved as a cycle.
    @Test
    public void testUnresolvedCycle()
    {
        PairwiseMargins margins = new Margins()
                .challenge("v1", "v2", 40).challenge("v2", "v3", 40).challenge("v3", "v1", 40).build();

        List<OncologyGroupRepresentativeSelection> selections =
                select(support(candidate("v1", 100), candidate("v2", 98), candidate("v3", 96)), margins);

        assertEquals(OncologyGroupOutcome.CYCLE, group(selections, GROUP_A).outcome());
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig decisively challenges the leader: unresolved, and the leader is left contested.
    @Test
    public void testUnresolvedMinorChallenger()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate("v1", 100), candidate("v3", 10)), new Margins().challenge("v3", "v1", 40).build());

        assertEquals(OncologyGroupOutcome.MINOR_RIVAL, group(selections, GROUP_A).outcome());
        assertEquals(ContigRole.CONTESTED, role(selections, "v1"));
        assertEquals(ContigRole.MINOR_CHALLENGER, role(selections, "v3"));
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig that challenges nobody: the leader is still chosen and the low-vote contig is minor.
    @Test
    public void testResolvedWithMinorBystander()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate("v1", 100), candidate("v3", 10)), noChallenges());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(ContigRole.MINOR, role(selections, "v3"));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
    }

    // A single candidate resolves trivially.
    @Test
    public void testSoleContig()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate("v1", 100)), noChallenges());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(OncologyGroupOutcome.ONE_CANDIDATE, group(selections, GROUP_A).outcome());
    }

    // A tenth of the group's reads winning by a decisive margin meets the challenge threshold of a tenth.
    // The denominator is every read aligning to the group, so here 200 votes over both contigs means 20 reads.
    @Test
    public void testChallengeAtThreshold()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate("v1", 100), candidate("v2", 100)), new Margins().challenge("v1", "v2", 20).build());

        assertEquals(ContigRole.SECONDARY, role(selections, "v2"));
    }

    @Test
    public void testNoChallengeBelowThreshold()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate("v1", 100), candidate("v2", 100)), new Margins().challenge("v1", "v2", 19).build());

        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, "v2"));
    }

    // Plenty of winning reads, but none of them wins by enough bases to count.
    @Test
    public void testNoChallengeBelowMargin()
    {
        PairwiseMargins margins = new Margins().challengeAtMargin("v1", "v2", MIN_CHALLENGE_MARGIN - 1, 100).build();

        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate("v1", 100), candidate("v2", 100)), margins);

        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, "v2"));
    }

    // Contigs the prefilter rejected are reported but take no part, and a group with none of them keeps nothing.
    @Test
    public void testRejectedContigsTakeNoPart()
    {
        List<ContigSupport> support = support(
                candidate("v1", 100),
                rejected("v2", ContigFilterStatus.LOW_VOTE_DENSITY),
                rejected("h1", ContigFilterStatus.LOW_COVERAGE));

        List<OncologyGroupRepresentativeSelection> selections = select(support, noChallenges());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, "v1"));
        assertEquals(List.of("v2"), names(group(selections, GROUP_A).rejected()));
        assertEquals(OncologyGroupOutcome.NO_CANDIDATES, group(selections, GROUP_H).outcome());
        assertEquals(OncologyGroupResolution.NO_CANDIDATES, group(selections, GROUP_H).resolution());
        assertNull(group(selections, GROUP_H).representative());
        assertTrue(group(selections, GROUP_H).candidates().isEmpty());
    }

    // Rank follows the votes order, whatever order the contigs arrive in.
    @Test
    public void testVotesRankFollowsSupport()
    {
        List<OncologyGroupRepresentativeSelection> selections =
                select(support(candidate("v2", 95), candidate("v1", 100)), noChallenges());

        assertEquals(1, candidate(selections, GROUP_A, "v1").votesRank());
        assertEquals(2, candidate(selections, GROUP_A, "v2").votesRank());
    }

    private static List<OncologyGroupRepresentativeSelection> select(List<ContigSupport> support, PairwiseMargins margins)
    {
        return RepresentativeSelector.select(support, margins, groupReadCounts(support));
    }

    // Each group is given as many reads as its contigs have votes, so a contig's votes read as its share of the group.
    private static Map<OncologyGroup, Integer> groupReadCounts(List<ContigSupport> support)
    {
        Map<OncologyGroup, Integer> readCounts = new HashMap<>();
        support.forEach(contig -> readCounts.merge(
                contig.contig().oncologyGroup(), (int) Math.round(contig.readVotes()), Integer::sum));
        return readCounts;
    }

    private static RepresentativeCandidate candidate(
            List<OncologyGroupRepresentativeSelection> selections, OncologyGroup oncologyGroup, String contig)
    {
        return group(selections, oncologyGroup).candidates().stream()
                .filter(candidate -> candidate.contig().equals(REFERENCE.contig(contig)))
                .findFirst()
                .orElseThrow();
    }

    private static OncologyGroupRepresentativeSelection group(List<OncologyGroupRepresentativeSelection> selections,
            OncologyGroup oncologyGroup)
    {
        return selections.stream()
                .filter(selection -> selection.oncologyGroup().equals(oncologyGroup))
                .findFirst()
                .orElseThrow();
    }

    private static RepresentativeCandidate result(List<OncologyGroupRepresentativeSelection> selections, String contig)
    {
        return selections.stream()
                .flatMap(selection -> selection.candidates().stream())
                .filter(candidate -> candidate.contig().name().equals(contig))
                .findFirst()
                .orElseThrow();
    }

    private static ContigRole role(List<OncologyGroupRepresentativeSelection> selections, String contig)
    {
        return result(selections, contig).role();
    }

    private static List<String> names(List<ContigSupport> support)
    {
        return support.stream().map(contig -> contig.contig().name()).toList();
    }

    private static List<ContigSupport> support(ContigSupport... support)
    {
        return List.of(support);
    }

    private static ContigSupport candidate(String contig, double votes)
    {
        return contigSupport(contig, ContigFilterStatus.CANDIDATE, votes);
    }

    private static ContigSupport rejected(String contig, ContigFilterStatus filterStatus)
    {
        return contigSupport(contig, filterStatus, 100);
    }

    private static ContigSupport contigSupport(String contig, ContigFilterStatus filterStatus, double votes)
    {
        SummaryStats dummy = SummaryStats.from(new int[] { 1 });
        return new ContigSupport(REFERENCE.contig(contig), filterStatus, 100, 0, dummy, 0, LENGTH / 2, dummy, dummy, votes);
    }

    private static ViralReference reference(String... contigNames)
    {
        List<ViralContig> contigs = new ArrayList<>();
        List<SAMSequenceRecord> records = new ArrayList<>();
        for(String name : contigNames)
        {
            OncologyGroup group = name.startsWith("h") ? GROUP_H : GROUP_A;
            contigs.add(new ViralContig(name, LENGTH, "Virus " + name, group));
            records.add(new SAMSequenceRecord(name, LENGTH));
        }
        return new ViralReference(contigs, new SAMSequenceDictionary(records));
    }

    private static PairwiseMargins noChallenges()
    {
        return new Margins().build();
    }

    // Builds pairwise margins with controlled challenge edges: each challenge places its reads at a margin above the
    // selection margin, so they count as decisive.
    private static class Margins
    {
        private final Map<ContigPair, NavigableMap<Integer, Integer>> mMargins = new HashMap<>();
        private final Map<ContigPair, Integer> mShared = new HashMap<>();

        private Margins challenge(String subject, String opponent, int reads)
        {
            return challengeAtMargin(subject, opponent, MIN_CHALLENGE_MARGIN * 2, reads);
        }

        private Margins challengeAtMargin(String subject, String opponent, int margin, int reads)
        {
            ContigPair pair = new ContigPair(REFERENCE.contig(subject), REFERENCE.contig(opponent));
            mMargins.computeIfAbsent(pair, key -> new TreeMap<>()).merge(margin, reads, Integer::sum);
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
