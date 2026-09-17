package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.emptyList;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.COMPARABLE_VOTE_RATIO;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_READS;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE_LOWER;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_VOTES_PER_BASE;

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
    public RepresentativeSelectionResult classify(
            Map<String, ContigStats> contigStats, PairwiseMargins pairwise, ViralReference reference)
    {
        Map<String, List<ContigStats>> contigsByOncologyGroup = contigStats.values().stream()
                .collect(groupingBy(stats -> reference.contig(stats.contig()).oncologyGroup()));

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
        Set<String> coveredContigs = contigNames(covered);
        List<ContigClassification> lowCoverage = groupContigs.stream()
                .filter(stats -> !coveredContigs.contains(stats.contig()))
                .map(stats -> droppedClassification(stats, oncologyGroup, ContigFilterStatus.LOW_COVERAGE))
                .toList();

        if(covered.isEmpty())
        {
            return new OncologyGroupResult(oncologyGroup, null, lowCoverage);
        }

        double voteTotal = covered.stream().mapToDouble(ContigStats::readVotes).sum();
        Map<String, Integer> votesRankByContig = votesRanks(covered);

        List<ContigStats> candidates = covered.stream()
                .filter(stats -> passesVoteDensity(stats, pairwise.meanReadLength()))
                .toList();
        Set<String> candidateContigs = contigNames(candidates);
        List<ContigClassification> lowVoteDensity = covered.stream()
                .filter(stats -> !candidateContigs.contains(stats.contig()))
                .map(stats -> lowVoteDensityClassification(stats, oncologyGroup, votesRankByContig, voteTotal))
                .toList();

        if(candidates.isEmpty())
        {
            List<ContigClassification> contigClassifications = Stream.concat(lowCoverage.stream(), lowVoteDensity.stream()).toList();
            return new OncologyGroupResult(oncologyGroup, voteTotal, contigClassifications);
        }

        ChallengeResolution resolution = resolveByChallenges(candidates, votesRankByContig, voteTotal, pairwise);
        double topVoteShare = candidates.stream().mapToDouble(stats -> voteShare(stats.readVotes(), voteTotal)).max().orElse(0.0);
        List<ContigClassification> candidateClassifications = candidates.stream()
                .map(candidate -> candidateClassification(
                        candidate, oncologyGroup, resolution, candidates, votesRankByContig, voteTotal, topVoteShare, pairwise))
                .toList();

        List<ContigClassification> contigClassifications =
                Stream.of(lowCoverage, lowVoteDensity, candidateClassifications).flatMap(List::stream).toList();

        return new OncologyGroupResult(oncologyGroup, voteTotal, contigClassifications);
    }

    // Builds the challenge graph over an oncology group's candidates and reduces it to per-contig roles and an outcome.
    private ChallengeResolution resolveByChallenges(
            List<ContigStats> candidates, Map<String, Integer> votesRankByContig, double voteTotal, PairwiseMargins pairwise)
    {
        List<String> contigs = contigNames(candidates).stream().toList();
        if(contigs.size() == 1)
        {
            return new ChallengeResolution(
                    Map.of(contigs.get(0), ContigRole.REPRESENTATIVE), Set.of(contigs.get(0)), OncologyGroupSubOutcome.ONE_CANDIDATE);
        }

        Set<String> comparable = abundantContigs(candidates, voteTotal);

        boolean minorChallengesAbundant = contigs.stream()
                .filter(contig -> !comparable.contains(contig))
                .anyMatch(minor -> comparable.stream().anyMatch(peer -> challenges(minor, peer, pairwise, voteTotal)));

        Map<String, Boolean> challengedByPeer = new HashMap<>();
        for(String contig : comparable)
        {
            challengedByPeer.put(
                    contig,
                    comparable.stream().anyMatch(peer -> !peer.equals(contig) && challenges(peer, contig, pairwise, voteTotal)));
        }
        List<String> unchallenged = comparable.stream().filter(contig -> !challengedByPeer.get(contig)).toList();

        String representative = null;
        OncologyGroupSubOutcome subOutcome;
        if(minorChallengesAbundant)
        {
            subOutcome = OncologyGroupSubOutcome.MINOR_RIVAL;
        }
        else if(unchallenged.isEmpty())
        {
            subOutcome = hasMutualChallenge(comparable, pairwise, voteTotal)
                    ? OncologyGroupSubOutcome.MUTUAL : OncologyGroupSubOutcome.CYCLE;
        }
        else
        {
            subOutcome = OncologyGroupSubOutcome.RESOLVED_CANDIDATES;
            representative = unchallenged.stream().min(Comparator.comparingInt(votesRankByContig::get)).orElseThrow();
        }

        Map<String, ContigRole> roles = new HashMap<>();
        for(String contig : contigs)
        {
            roles.put(contig, roleFor(contig, representative, comparable, challengedByPeer, pairwise, voteTotal));
        }
        return new ChallengeResolution(roles, comparable, subOutcome);
    }

    private ContigRole roleFor(
            String contig, String representative, Set<String> comparable,
            Map<String, Boolean> challengedByPeer, PairwiseMargins pairwise, double voteTotal)
    {
        if(contig.equals(representative))
        {
            return ContigRole.REPRESENTATIVE;
        }
        if(!comparable.contains(contig))
        {
            boolean challengesAbundant = comparable.stream().anyMatch(peer -> challenges(contig, peer, pairwise, voteTotal));
            return challengesAbundant ? ContigRole.MINOR_CHALLENGER : ContigRole.MINOR;
        }
        if(challengedByPeer.get(contig))
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
        return stats.readVotes() >= voteFloorPerBase * stats.contigLength();
    }

    // Contigs whose vote share is near the oncology group's top: the abundant contenders that contest each other.
    private static Set<String> abundantContigs(List<ContigStats> candidates, double voteTotal)
    {
        double topVoteShare = candidates.stream().mapToDouble(stats -> voteShare(stats.readVotes(), voteTotal)).max().orElse(0.0);
        return candidates.stream()
                .filter(stats -> topVoteShare > 0 && voteShare(stats.readVotes(), voteTotal) >= COMPARABLE_VOTE_RATIO * topVoteShare)
                .map(ContigStats::contig)
                .collect(toSet());
    }

    private static boolean hasMutualChallenge(Set<String> comparable, PairwiseMargins pairwise, double voteTotal)
    {
        return comparable.stream().anyMatch(subject -> comparable.stream().anyMatch(opponent -> !subject.equals(opponent)
                && challenges(subject, opponent, pairwise, voteTotal) && challenges(opponent, subject, pairwise, voteTotal)));
    }

    // A contig challenges another when a high enough fraction of the oncology group's reads prefer it by at least the margin.
    private static boolean challenges(String subject, String opponent, PairwiseMargins pairwise, double voteTotal)
    {
        if(voteTotal <= 0)
        {
            return false;
        }
        return pairwise.challengeReads(subject, opponent, MIN_CHALLENGE_MARGIN) / voteTotal >= MIN_CHALLENGE_READS;
    }

    private ContigClassification candidateClassification(
            ContigStats candidate, String oncologyGroup, ChallengeResolution resolution, List<ContigStats> candidates,
            Map<String, Integer> votesRankByContig, double voteTotal, double topVoteShare, PairwiseMargins pairwise)
    {
        double candidateVoteShare = voteShare(candidate.readVotes(), voteTotal);
        List<Integer> challenges = candidateRanks(
                candidates, candidate, other -> challenges(candidate.contig(), other, pairwise, voteTotal), votesRankByContig);
        List<Integer> challengedBy = candidateRanks(
                candidates, candidate, other -> challenges(other, candidate.contig(), pairwise, voteTotal), votesRankByContig);

        return new ContigClassification(
                candidate.contig(), oncologyGroup, ContigFilterStatus.CANDIDATE,
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
            Predicate<String> matches, Map<String, Integer> votesRankByContig)
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
            ContigStats stats, String oncologyGroup, Map<String, Integer> votesRankByContig, double voteTotal)
    {
        return new ContigClassification(
                stats.contig(), oncologyGroup, ContigFilterStatus.LOW_VOTE_DENSITY,
                votesRankByContig.get(stats.contig()), voteShare(stats.readVotes(), voteTotal),
                null, null, emptyList(), emptyList(), null, null, null);
    }

    private static ContigClassification droppedClassification(ContigStats stats, String oncologyGroup, ContigFilterStatus status)
    {
        return new ContigClassification(
                stats.contig(), oncologyGroup, status, null, null, null, null, emptyList(), emptyList(), null, null, null);
    }

    // Rank 1 = most read votes; contig name breaks ties for determinism.
    private static Map<String, Integer> votesRanks(List<ContigStats> covered)
    {
        List<ContigStats> ordered = covered.stream()
                .sorted(Comparator.comparingDouble(ContigStats::readVotes).reversed().thenComparing(ContigStats::contig))
                .toList();
        Map<String, Integer> votesRankByContig = new HashMap<>();
        for(int i = 0; i < ordered.size(); ++i)
        {
            votesRankByContig.put(ordered.get(i).contig(), i + 1);
        }
        return votesRankByContig;
    }

    private static Set<String> contigNames(List<ContigStats> stats)
    {
        return stats.stream().map(ContigStats::contig).collect(toSet());
    }

    private static double voteShare(double votes, double voteTotal)
    {
        return voteTotal > 0 ? votes / voteTotal : 0.0;
    }

    private record ChallengeResolution(
            Map<String, ContigRole> roles,
            Set<String> comparable,
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
