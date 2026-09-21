package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.junit.Test;

public class RepresentativeSelectorTest
{
    private static final int LENGTH = 1000;
    private static final double MEAN_READ_LENGTH = 150.0;

    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_H = new OncologyGroup("Group H");

    private static final ViralContig V1 = new ViralContig("v1", LENGTH, "Virus v1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", LENGTH, "Virus v2", GROUP_A);
    private static final ViralContig V3 = new ViralContig("v3", LENGTH, "Virus v3", GROUP_A);
    private static final ViralContig H1 = new ViralContig("h1", LENGTH, "Virus h1", GROUP_H);

    // Two contigs with similar votes and no challenge between them: the one with more votes leads, the other is its twin.
    @Test
    public void testResolvedTwins()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate(V1, 100), candidate(V2, 95)), margins());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, V1));
        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, V2));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
        assertEquals(OncologyGroupResolution.RESOLVED, group(selections, GROUP_A).resolution());
        assertEquals(V1, group(selections, GROUP_A).representative());
    }

    // The leader decisively challenges the other contig: the other becomes secondary and the leader is still chosen.
    @Test
    public void testResolvedWithSecondary()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 95)), margins(challenge(V1, V2, 40)));

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, V1));
        assertEquals(ContigRole.SECONDARY, role(selections, V2));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
        assertEquals(Set.of(V2), result(selections, V1).challenges());
        assertEquals(Set.of(V1), result(selections, V2).challengedBy());
    }

    // Two contigs challenge each other, so neither can be chosen: unresolved, no representative.
    @Test
    public void testUnresolvedMutual()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 95)), margins(challenge(V1, V2, 40), challenge(V2, V1, 40)));

        assertEquals(OncologyGroupOutcome.MUTUAL, group(selections, GROUP_A).outcome());
        assertEquals(OncologyGroupResolution.UNRESOLVED, group(selections, GROUP_A).resolution());
        assertNull(group(selections, GROUP_A).representative());
    }

    // Three contigs challenge in a loop (v1 beats v2 beats v3 beats v1), so none is unchallenged: unresolved as a cycle.
    @Test
    public void testUnresolvedCycle()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 98), candidate(V3, 96)),
                margins(challenge(V1, V2, 40), challenge(V2, V3, 40), challenge(V3, V1, 40)));

        assertEquals(OncologyGroupOutcome.CYCLE, group(selections, GROUP_A).outcome());
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig decisively challenges the leader: unresolved, and the leader is left contested.
    @Test
    public void testUnresolvedMinorChallenger()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V3, 10)), margins(challenge(V3, V1, 40)));

        assertEquals(OncologyGroupOutcome.MINOR_RIVAL, group(selections, GROUP_A).outcome());
        assertEquals(ContigRole.CONTESTED, role(selections, V1));
        assertEquals(ContigRole.MINOR_CHALLENGER, role(selections, V3));
        assertNull(group(selections, GROUP_A).representative());
    }

    // A low-vote contig that challenges nobody: the leader is still chosen and the low-vote contig is minor.
    @Test
    public void testResolvedWithMinorBystander()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate(V1, 100), candidate(V3, 10)), margins());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, V1));
        assertEquals(ContigRole.MINOR, role(selections, V3));
        assertEquals(OncologyGroupOutcome.RESOLVED_CANDIDATES, group(selections, GROUP_A).outcome());
    }

    // A single candidate resolves trivially.
    @Test
    public void testSoleContig()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate(V1, 100)), margins());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, V1));
        assertEquals(OncologyGroupOutcome.ONE_CANDIDATE, group(selections, GROUP_A).outcome());
    }

    // A tenth of the group's reads winning by a decisive margin meets the challenge threshold of a tenth.
    @Test
    public void testChallengeAtThreshold()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 100)), margins(challenge(V1, V2, 20)), 200);

        assertEquals(ContigRole.SECONDARY, role(selections, V2));
    }

    @Test
    public void testNoChallengeBelowThreshold()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 100)), margins(challenge(V1, V2, 19)), 200);

        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, V2));
    }

    // Plenty of winning reads, but none of them wins by enough bases to count.
    @Test
    public void testNoChallengeBelowMargin()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(
                support(candidate(V1, 100), candidate(V2, 100)),
                margins(challenge(V1, V2, 100, MIN_CHALLENGE_MARGIN - 1)));

        assertEquals(ContigRole.REPRESENTATIVE_TWIN, role(selections, V2));
    }

    // Contigs the prefilter rejected are reported but take no part in selection, and a group with none of them keeps nothing.
    @Test
    public void testRejectedContigsExcluded()
    {
        List<ContigSupport> support = support(
                candidate(V1, 100),
                rejected(V2, ContigFilterStatus.LOW_VOTE_DENSITY),
                rejected(H1, ContigFilterStatus.LOW_COVERAGE));

        List<OncologyGroupRepresentativeSelection> selections = select(support, margins());

        assertEquals(ContigRole.REPRESENTATIVE, role(selections, V1));
        assertEquals(List.of(V2), contigs(group(selections, GROUP_A).rejected()));
        assertEquals(OncologyGroupOutcome.NO_CANDIDATES, group(selections, GROUP_H).outcome());
        assertEquals(OncologyGroupResolution.NO_CANDIDATES, group(selections, GROUP_H).resolution());
        assertNull(group(selections, GROUP_H).representative());
        assertTrue(group(selections, GROUP_H).candidates().isEmpty());
    }

    // Rank follows the votes order, whatever order the contigs arrive in.
    @Test
    public void testVotesRankFollowsSupport()
    {
        List<OncologyGroupRepresentativeSelection> selections = select(support(candidate(V2, 95), candidate(V1, 100)), margins());

        assertEquals(1, result(selections, V1).votesRank());
        assertEquals(2, result(selections, V2).votesRank());
    }

    // Group read count is derived from sum of votes, which is plausible enough.
    // Use the other overload if you need to set the read count precisely.
    private static List<OncologyGroupRepresentativeSelection> select(List<ContigSupport> support, PairwiseMargins margins)
    {
        Map<OncologyGroup, Integer> readCounts = new HashMap<>();
        support.forEach(contig -> readCounts.merge(
                contig.contig().oncologyGroup(), (int) Math.round(contig.readVotes()), Integer::sum));
        return RepresentativeSelector.select(support, margins, readCounts);
    }

    private static List<OncologyGroupRepresentativeSelection> select(
            List<ContigSupport> support, PairwiseMargins margins, int groupReads)
    {
        Map<OncologyGroup, Integer> readCounts = support.stream()
                .map(contig -> contig.contig().oncologyGroup())
                .distinct()
                .collect(toMap(group -> group, group -> groupReads));
        return RepresentativeSelector.select(support, margins, readCounts);
    }

    @SafeVarargs
    private static PairwiseMargins margins(List<ViralAlignment>... challenges)
    {
        List<ViralAlignment> alignments = Stream.of(challenges).flatMap(List::stream).toList();
        return PairwiseMargins.from(ViralAlignments.from(alignments, MEAN_READ_LENGTH));
    }

    private static List<ViralAlignment> challenge(ViralContig subject, ViralContig opponent, int reads)
    {
        return challenge(subject, opponent, reads, MIN_CHALLENGE_MARGIN * 2);
    }

    // Reads the subject fits `margin` bases better than the opponent.
    // Read names are unique per pair, so separate challenges share no read.
    private static List<ViralAlignment> challenge(ViralContig subject, ViralContig opponent, int reads, int margin)
    {
        List<ViralAlignment> alignments = new ArrayList<>();
        for(int i = 0; i < reads; ++i)
        {
            String readName = subject.name() + opponent.name() + i;
            alignments.add(alignment(readName, subject, 0));
            alignments.add(alignment(readName, opponent, margin));
        }
        return alignments;
    }

    private static ViralAlignment alignment(String readName, ViralContig contig, int divergence)
    {
        return new ViralAlignment(
                readName, contig, 1, LENGTH, 0, 0, 100, divergence, List.of(new AlignedInterval(1, LENGTH)));
    }

    private static OncologyGroupRepresentativeSelection group(List<OncologyGroupRepresentativeSelection> selections,
            OncologyGroup oncologyGroup)
    {
        return selections.stream()
                .filter(selection -> selection.oncologyGroup().equals(oncologyGroup))
                .findFirst()
                .orElseThrow();
    }

    private static RepresentativeCandidate result(List<OncologyGroupRepresentativeSelection> selections, ViralContig contig)
    {
        return selections.stream()
                .flatMap(selection -> selection.candidates().stream())
                .filter(candidate -> candidate.contig().equals(contig))
                .findFirst()
                .orElseThrow();
    }

    private static ContigRole role(List<OncologyGroupRepresentativeSelection> selections, ViralContig contig)
    {
        return result(selections, contig).role();
    }

    private static List<ViralContig> contigs(List<ContigSupport> support)
    {
        return support.stream().map(ContigSupport::contig).toList();
    }

    private static List<ContigSupport> support(ContigSupport... support)
    {
        return List.of(support);
    }

    private static ContigSupport candidate(ViralContig contig, double votes)
    {
        return contigSupport(contig, ContigFilterStatus.CANDIDATE, votes);
    }

    private static ContigSupport rejected(ViralContig contig, ContigFilterStatus filterStatus)
    {
        return contigSupport(contig, filterStatus, 100);
    }

    // Selection reads only the contig, its filter status and its read votes. The remaining fields are filler, set to
    // values consistent with a present contig so nothing reads as contradictory.
    private static ContigSupport contigSupport(ViralContig contig, ContigFilterStatus filterStatus, double votes)
    {
        SummaryStats filler = SummaryStats.from(new int[] { 1 });
        return new ContigSupport(contig, filterStatus, 100, 0, filler, 0, LENGTH / 2, filler, filler, votes);
    }
}
