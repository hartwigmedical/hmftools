package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.emptyList;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;
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
// Then a pairwise "challenges" graph is generated to determine the presence of "rival" contigs - contigs with low
// overall support, but are decisively supported by a subset of reads. I.e. 1 contig doesn't explain the whole viral
// genome in the sample.
public class RepresentativeSelector
{
    public RepresentativeSelectionResult classify(Collection<ContigStats> contigStats, PairwiseMargins pairwise)
    {
        Map<String, List<ContigStats>> contigsByOncologyGroup = contigStats.stream()
                .collect(groupingBy(stats -> stats.contig().oncologyGroup()));

        List<OncologyGroupResult> groupResults = contigsByOncologyGroup.entrySet().stream()
                .map(entry -> classifyOncologyGroup(entry.getKey(), entry.getValue(), pairwise))
                .toList();

        List<ContigClassification> classifications = groupResults.stream()
                .flatMap(result -> result.classifications().stream())
                .toList();
        Map<String, Double> oncologyGroupVoteTotals = groupResults.stream()
                .filter(result -> result.voteTotal() != null)
                .collect(toMap(OncologyGroupResult::oncologyGroup, OncologyGroupResult::voteTotal));

        return new RepresentativeSelectionResult(classifications, oncologyGroupVoteTotals);
    }

    private OncologyGroupResult classifyOncologyGroup(String oncologyGroup, List<ContigStats> groupContigs, PairwiseMargins pairwise)
    {
        List<ContigStats> covered = contigsPassingCoverage(groupContigs);
        Set<ViralContig> coveredContigs = contigSet(covered);
        List<ContigClassification> lowCoverage = groupContigs.stream()
                .filter(stats -> !coveredContigs.contains(stats.contig()))
                .map(RepresentativeSelector::lowCoverageClassification)
                .toList();

        if(covered.isEmpty())
        {
            return new OncologyGroupResult(oncologyGroup, null, lowCoverage);
        }

        double voteTotal = covered.stream().mapToDouble(ContigStats::readVotes).sum();
        Map<ViralContig, Integer> votesRankByContig = votesRanks(covered);

        List<ContigStats> candidates = covered.stream()
                .filter(stats -> passesVoteDensity(stats, pairwise.meanReadLength()))
                .toList();
        Set<ViralContig> candidateContigs = contigSet(candidates);
        List<ContigClassification> lowVoteDensity = covered.stream()
                .filter(stats -> !candidateContigs.contains(stats.contig()))
                .map(stats -> lowVoteDensityClassification(stats, votesRankByContig, voteTotal))
                .toList();

        if(candidates.isEmpty())
        {
            List<ContigClassification> contigClassifications = Stream.concat(lowCoverage.stream(), lowVoteDensity.stream()).toList();
            return new OncologyGroupResult(oncologyGroup, voteTotal, contigClassifications);
        }

        ChallengeGraph graph = ChallengeGraph.build(candidates, pairwise, voteTotal);
        ChallengeResolution resolution = resolveByChallenges(votesRankByContig, graph);
        List<ContigClassification> candidateClassifications = candidates.stream()
                .map(candidate -> candidateClassification(candidate, resolution, candidates, votesRankByContig, voteTotal, graph))
                .toList();

        List<ContigClassification> contigClassifications =
                Stream.of(lowCoverage, lowVoteDensity, candidateClassifications).flatMap(List::stream).toList();

        return new OncologyGroupResult(oncologyGroup, voteTotal, contigClassifications);
    }

    // Reduces the challenge graph over an oncology group's candidates to per-contig roles and an outcome.
    private ChallengeResolution resolveByChallenges(Map<ViralContig, Integer> votesRankByContig, ChallengeGraph graph)
    {
        List<ViralContig> contigs = graph.contigs();
        if(contigs.size() == 1)
        {
            return new ChallengeResolution(
                    Map.of(contigs.get(0), ContigRole.REPRESENTATIVE), Set.of(contigs.get(0)), OncologyGroupSubOutcome.ONE_CANDIDATE);
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
        OncologyGroupSubOutcome subOutcome;
        if(minorChallengesAbundant)
        {
            subOutcome = OncologyGroupSubOutcome.MINOR_RIVAL;
        }
        else if(unchallenged.isEmpty())
        {
            subOutcome = hasMutualChallenge(comparable, graph) ? OncologyGroupSubOutcome.MUTUAL : OncologyGroupSubOutcome.CYCLE;
        }
        else
        {
            subOutcome = OncologyGroupSubOutcome.RESOLVED_CANDIDATES;
            representative = unchallenged.stream().min(Comparator.comparingInt(votesRankByContig::get)).orElseThrow();
        }

        Map<ViralContig, ContigRole> roles = new HashMap<>();
        for(ViralContig contig : contigs)
        {
            roles.put(contig, roleFor(contig, representative, comparable, challengedByPeer, graph));
        }
        return new ChallengeResolution(roles, comparable, subOutcome);
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
            ContigStats candidate, ChallengeResolution resolution, List<ContigStats> candidates,
            Map<ViralContig, Integer> votesRankByContig, double voteTotal, ChallengeGraph graph)
    {
        double candidateVoteShare = voteShare(candidate.readVotes(), voteTotal);
        double topVoteShare = graph.topVoteShare();
        List<Integer> challenges = candidateRanks(
                candidates, candidate, other -> graph.challenges(candidate.contig(), other), votesRankByContig);
        List<Integer> challengedBy = candidateRanks(
                candidates, candidate, other -> graph.challenges(other, candidate.contig()), votesRankByContig);

        return new ContigClassification(
                candidate.contig(), ContigFilterStatus.CANDIDATE,
                votesRankByContig.get(candidate.contig()), candidateVoteShare,
                topVoteShare > 0 ? candidateVoteShare / topVoteShare : 0.0,
                resolution.comparable().contains(candidate.contig()),
                challenges, challengedBy,
                resolution.roles().get(candidate.contig()),
                resolution.subOutcome().outcome(), resolution.subOutcome());
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
            ContigStats stats, Map<ViralContig, Integer> votesRankByContig, double voteTotal)
    {
        return new ContigClassification(
                stats.contig(), ContigFilterStatus.LOW_VOTE_DENSITY,
                votesRankByContig.get(stats.contig()), voteShare(stats.readVotes(), voteTotal),
                null, null, emptyList(), emptyList(), null, null, null);
    }

    private static ContigClassification lowCoverageClassification(ContigStats stats)
    {
        return new ContigClassification(
                stats.contig(), ContigFilterStatus.LOW_COVERAGE, null, null, null, null, emptyList(), emptyList(), null, null, null);
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

    private static double voteShare(double votes, double voteTotal)
    {
        return voteTotal > 0 ? votes / voteTotal : 0.0;
    }

    private record ChallengeResolution(
            Map<ViralContig, ContigRole> roles,
            Set<ViralContig> comparable,
            OncologyGroupSubOutcome subOutcome
    )
    {
    }

    private record OncologyGroupResult(
            String oncologyGroup,
            @Nullable Double voteTotal,
            List<ContigClassification> classifications
    )
    {
    }
}
