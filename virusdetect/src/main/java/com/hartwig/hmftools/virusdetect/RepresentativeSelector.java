package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.emptyList;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE_LOWER;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_VOTES_PER_BASE;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

// Per oncology group, pick at most one representative contig.
// First, contigs are filtered on coverage and read votes.
// Then a "challenges" graph determines the presence of "rival" contigs - contigs with low overall support, but
// decisively supported by a subset of reads. I.e. 1 contig doesn't explain the whole viral genome in the sample.
public class RepresentativeSelector
{
    public RepresentativeSelectionResult classify(
            Collection<ContigStats> contigStats, PairwiseMargins margins, double meanReadLength)
    {
        if(meanReadLength <= 0)
        {
            throw new IllegalArgumentException("invalid mean read length: " + meanReadLength);
        }

        List<ContigClassification> classifications = contigStats.stream()
                .collect(groupingBy(stats -> stats.contig().oncologyGroup()))
                .values().stream()
                .flatMap(groupContigs -> classifyOncologyGroup(groupContigs, margins, meanReadLength).stream())
                .toList();

        return new RepresentativeSelectionResult(classifications);
    }

    private List<ContigClassification> classifyOncologyGroup(
            List<ContigStats> groupContigs, PairwiseMargins margins, double meanReadLength)
    {
        double allVotesTotal = groupContigs.stream().mapToDouble(ContigStats::readVotes).sum();

        List<ContigStats> covered = contigsPassingCoverage(groupContigs);
        Set<ViralContig> coveredContigs = contigSet(covered);

        if(covered.isEmpty())
        {
            return groupContigs.stream().map(stats -> lowCoverageClassification(stats, allVotesTotal, 0.0)).toList();
        }

        Map<ViralContig, Integer> votesRankByContig = votesRanks(covered);

        List<ContigStats> candidates = covered.stream()
                .filter(stats -> passesVoteDensity(stats, meanReadLength))
                .toList();
        Set<ViralContig> candidateContigs = contigSet(candidates);
        double candidateVotesTotal = candidates.stream().mapToDouble(ContigStats::readVotes).sum();
        double topCandidateVotes = candidates.stream().mapToDouble(ContigStats::readVotes).max().orElse(0.0);

        List<ContigClassification> lowCoverage = groupContigs.stream()
                .filter(stats -> !coveredContigs.contains(stats.contig()))
                .map(stats -> lowCoverageClassification(stats, allVotesTotal, candidateVotesTotal))
                .toList();
        List<ContigClassification> lowVoteDensity = covered.stream()
                .filter(stats -> !candidateContigs.contains(stats.contig()))
                .map(stats -> lowVoteDensityClassification(stats, votesRankByContig, allVotesTotal, candidateVotesTotal))
                .toList();

        if(candidates.isEmpty())
        {
            return Stream.concat(lowCoverage.stream(), lowVoteDensity.stream()).toList();
        }

        ChallengeGraph graph = ChallengeGraph.build(candidates, groupContigs, margins);
        ChallengeResolution challengeResolution = resolveByChallenges(votesRankByContig, graph);
        List<ContigClassification> candidateClassifications = candidates.stream()
                .map(candidate -> candidateClassification(
                        candidate, challengeResolution, candidates, votesRankByContig, allVotesTotal, candidateVotesTotal,
                        topCandidateVotes, graph))
                .toList();

        return Stream.of(lowCoverage, lowVoteDensity, candidateClassifications).flatMap(List::stream).toList();
    }

    // Reduces the challenge graph over an oncology group's candidates to per-contig roles and an outcome.
    private ChallengeResolution resolveByChallenges(Map<ViralContig, Integer> votesRankByContig, ChallengeGraph graph)
    {
        List<ViralContig> contigs = graph.contigs();
        if(contigs.size() == 1)
        {
            return new ChallengeResolution(
                    Map.of(contigs.get(0), ContigRole.REPRESENTATIVE), Set.of(contigs.get(0)), OncologyGroupOutcome.ONE_CANDIDATE);
        }

        Set<ViralContig> comparable = graph.comparable();

        boolean minorChallengesAbundant = contigs.stream()
                .filter(contig -> !comparable.contains(contig))
                .anyMatch(minor -> comparable.stream().anyMatch(peer -> graph.challenges(minor, peer)));

        Set<ViralContig> challengedByPeer = comparable.stream()
                .filter(contig -> comparable.stream().anyMatch(peer -> !peer.equals(contig) && graph.challenges(peer, contig)))
                .collect(toSet());
        List<ViralContig> unchallenged = comparable.stream().filter(contig -> !challengedByPeer.contains(contig)).toList();

        ViralContig representative = null;
        OncologyGroupOutcome outcome;
        if(minorChallengesAbundant)
        {
            outcome = OncologyGroupOutcome.MINOR_RIVAL;
        }
        else if(unchallenged.isEmpty())
        {
            outcome = hasMutualChallenge(comparable, graph) ? OncologyGroupOutcome.MUTUAL : OncologyGroupOutcome.CYCLE;
        }
        else
        {
            outcome = OncologyGroupOutcome.RESOLVED_CANDIDATES;
            representative = unchallenged.stream().min(Comparator.comparingInt(votesRankByContig::get)).orElseThrow();
        }

        Map<ViralContig, ContigRole> roles = new HashMap<>();
        for(ViralContig contig : contigs)
        {
            roles.put(contig, roleFor(contig, representative, comparable, challengedByPeer, graph));
        }
        return new ChallengeResolution(roles, comparable, outcome);
    }

    private ContigRole roleFor(
            ViralContig contig, @Nullable ViralContig representative, Set<ViralContig> comparable,
            Set<ViralContig> challengedByPeer, ChallengeGraph graph)
    {
        if(contig.equals(representative))
        {
            return ContigRole.REPRESENTATIVE;
        }
        if(!comparable.contains(contig))
        {
            boolean challengesAbundant = comparable.stream().anyMatch(peer -> graph.challenges(contig, peer));
            return challengesAbundant ? ContigRole.MINOR_CHALLENGER : ContigRole.MINOR;
        }
        if(challengedByPeer.contains(contig))
        {
            return ContigRole.SECONDARY;
        }
        // Abundant and unchallenged by peers, but not crowned: a resolved oncology group makes it a twin, else contested.
        return representative != null ? ContigRole.REPRESENTATIVE_TWIN : ContigRole.CONTESTED;
    }

    private static List<ContigStats> contigsPassingCoverage(List<ContigStats> groupContigs)
    {
        boolean present = groupContigs.stream().anyMatch(stats -> stats.coverageFraction() >= MIN_COVERAGE);
        if(!present)
        {
            return emptyList();
        }
        // If any contig is present, lower the coverage threshold to include contigs which straddle the threshold.
        return groupContigs.stream().filter(stats -> stats.coverageFraction() >= MIN_COVERAGE_LOWER).toList();
    }

    // Drops a contig that has enough coverage but few votes, because its reads align better to a different contig.
    private static boolean passesVoteDensity(ContigStats stats, double meanReadLength)
    {
        if(meanReadLength <= 0)
        {
            return true;
        }
        double voteFloorPerBase = MIN_VOTES_PER_BASE * MIN_COVERAGE / meanReadLength;
        return stats.readVotes() >= voteFloorPerBase * stats.contig().length();
    }

    private static boolean hasMutualChallenge(Set<ViralContig> comparable, ChallengeGraph graph)
    {
        return comparable.stream().anyMatch(subject -> comparable.stream().anyMatch(opponent ->
                !subject.equals(opponent) && graph.challenges(subject, opponent) && graph.challenges(opponent, subject)));
    }

    private ContigClassification candidateClassification(
            ContigStats candidate, ChallengeResolution challengeResolution, List<ContigStats> candidates,
            Map<ViralContig, Integer> votesRankByContig, double allVotesTotal, double candidateVotesTotal, double topCandidateVotes,
            ChallengeGraph graph)
    {
        List<Integer> challenges = candidateRanks(
                candidates, candidate, other -> graph.challenges(candidate.contig(), other), votesRankByContig);
        List<Integer> challengedBy = candidateRanks(
                candidates, candidate, other -> graph.challenges(other, candidate.contig()), votesRankByContig);

        return new ContigClassification(
                candidate.contig(), ContigFilterStatus.CANDIDATE,
                votesRankByContig.get(candidate.contig()),
                shareOrNull(candidate.readVotes(), allVotesTotal), shareOrNull(candidate.readVotes(), candidateVotesTotal),
                shareOrNull(candidate.readVotes(), topCandidateVotes),
                challengeResolution.comparable().contains(candidate.contig()),
                challenges, challengedBy,
                challengeResolution.roles().get(candidate.contig()),
                challengeResolution.outcome().resolution(), challengeResolution.outcome());
    }

    // The votes-ranks of the other candidates matching the challenge relation, sorted for stable output.
    private static List<Integer> candidateRanks(
            List<ContigStats> candidates, ContigStats subject,
            Predicate<ViralContig> matches, Map<ViralContig, Integer> votesRankByContig)
    {
        return candidates.stream()
                .map(ContigStats::contig)
                .filter(contig -> !contig.equals(subject.contig()))
                .filter(matches)
                .map(votesRankByContig::get)
                .sorted()
                .toList();
    }

    private static ContigClassification lowVoteDensityClassification(
            ContigStats stats, Map<ViralContig, Integer> votesRankByContig, double allVotesTotal, double candidateVotesTotal)
    {
        return new ContigClassification(
                stats.contig(), ContigFilterStatus.LOW_VOTE_DENSITY, votesRankByContig.get(stats.contig()),
                shareOrNull(stats.readVotes(), allVotesTotal), shareOrNull(stats.readVotes(), candidateVotesTotal),
                null, null, emptyList(), emptyList(), null, null, null);
    }

    private static ContigClassification lowCoverageClassification(ContigStats stats, double allVotesTotal, double candidateVotesTotal)
    {
        return new ContigClassification(
                stats.contig(), ContigFilterStatus.LOW_COVERAGE, null,
                shareOrNull(stats.readVotes(), allVotesTotal), shareOrNull(stats.readVotes(), candidateVotesTotal),
                null, null, emptyList(), emptyList(), null, null, null);
    }

    // Rank 1 = most read votes; contig name breaks ties for determinism.
    private static Map<ViralContig, Integer> votesRanks(List<ContigStats> covered)
    {
        List<ContigStats> ordered = covered.stream()
                .sorted(Comparator.comparingDouble(ContigStats::readVotes).reversed()
                        .thenComparing(stats -> stats.contig().name()))
                .toList();
        Map<ViralContig, Integer> votesRankByContig = new HashMap<>();
        for(int i = 0; i < ordered.size(); ++i)
        {
            votesRankByContig.put(ordered.get(i).contig(), i + 1);
        }
        return votesRankByContig;
    }

    private static Set<ViralContig> contigSet(List<ContigStats> stats)
    {
        return stats.stream().map(ContigStats::contig).collect(toSet());
    }

    @Nullable
    private static Double shareOrNull(double votes, double total)
    {
        return total > 0 ? votes / total : null;
    }

    private record ChallengeResolution(
            Map<ViralContig, ContigRole> roles,
            Set<ViralContig> comparable,
            OncologyGroupOutcome outcome
    )
    {
    }
}
