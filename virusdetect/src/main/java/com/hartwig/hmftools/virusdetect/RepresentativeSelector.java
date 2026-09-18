package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.emptyList;
import static java.util.stream.Collectors.groupingBy;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
        GroupCandidates prefiltered = GroupCandidates.from(groupContigs, meanReadLength);
        List<ContigStats> candidates = prefiltered.candidates();

        double allVotesTotal = groupContigs.stream().mapToDouble(ContigStats::readVotes).sum();
        double candidateVotesTotal = candidates.stream().mapToDouble(ContigStats::readVotes).sum();
        double topCandidateVotes = candidates.stream().mapToDouble(ContigStats::readVotes).max().orElse(0.0);

        List<ContigClassification> rejected = prefiltered.rejected().stream()
                .map(rejection -> rejectedClassification(rejection, allVotesTotal, candidateVotesTotal))
                .toList();

        if(candidates.isEmpty())
        {
            return rejected;
        }

        ChallengeGraph graph = ChallengeGraph.build(candidates, groupContigs, margins);
        RepresentativeChoice choice = RepresentativeChoice.from(candidates, graph);
        Map<ViralContig, Integer> votesRankByContig = votesRanks(candidates);

        List<ContigClassification> candidateClassifications = candidates.stream()
                .map(candidate -> candidateClassification(
                        candidate, choice, graph, candidates, votesRankByContig, allVotesTotal, candidateVotesTotal,
                        topCandidateVotes))
                .toList();

        return Stream.concat(rejected.stream(), candidateClassifications.stream()).toList();
    }

    private static ContigClassification candidateClassification(
            ContigStats candidate, RepresentativeChoice choice, ChallengeGraph graph, List<ContigStats> candidates,
            Map<ViralContig, Integer> votesRankByContig, double allVotesTotal, double candidateVotesTotal,
            double topCandidateVotes)
    {
        List<Integer> challenges = candidateRanks(
                candidates, candidate, other -> graph.challenges(candidate.contig(), other), votesRankByContig);
        List<Integer> challengedBy = candidateRanks(
                candidates, candidate, other -> graph.challenges(other, candidate.contig()), votesRankByContig);

        return new ContigClassification(
                candidate.contig(), ContigFilterStatus.CANDIDATE,
                votesRankByContig.get(candidate.contig()),
                shareOrNull(candidate.readVotes(), allVotesTotal),
                shareOrNull(candidate.readVotes(), candidateVotesTotal),
                shareOrNull(candidate.readVotes(), topCandidateVotes),
                graph.comparable().contains(candidate.contig()),
                challenges, challengedBy,
                choice.role(candidate.contig()),
                choice.outcome().resolution(), choice.outcome());
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

    private static ContigClassification rejectedClassification(
            GroupCandidates.Rejected rejection, double allVotesTotal, double candidateVotesTotal)
    {
        double votes = rejection.stats().readVotes();
        return new ContigClassification(
                rejection.stats().contig(), rejection.reason(), null,
                shareOrNull(votes, allVotesTotal), shareOrNull(votes, candidateVotesTotal),
                null, null, emptyList(), emptyList(), null, null, null);
    }

    // Rank 1 = best supported among the candidates.
    private static Map<ViralContig, Integer> votesRanks(List<ContigStats> candidates)
    {
        List<ContigStats> ordered = candidates.stream().sorted(ContigStats.BY_SUPPORT).toList();
        Map<ViralContig, Integer> votesRankByContig = new HashMap<>();
        for(int i = 0; i < ordered.size(); ++i)
        {
            votesRankByContig.put(ordered.get(i).contig(), i + 1);
        }
        return votesRankByContig;
    }

    @Nullable
    private static Double shareOrNull(double votes, double total)
    {
        return total > 0 ? votes / total : null;
    }
}
