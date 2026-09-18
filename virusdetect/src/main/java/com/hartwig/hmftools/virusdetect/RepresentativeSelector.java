package com.hartwig.hmftools.virusdetect;

import static java.util.Objects.requireNonNull;
import static java.util.stream.Collectors.groupingBy;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Stream;

// Per oncology group, pick at most one representative contig.
// First, contigs are filtered on coverage and read votes.
// Then a "challenges" graph determines the presence of "rival" contigs - contigs with low overall support, but
// decisively supported by a subset of reads. I.e. 1 contig doesn't explain the whole viral genome in the sample.
public class RepresentativeSelector
{
    public List<OncologyGroupSelection> select(
            Collection<ContigStats> contigStats, PairwiseMargins margins, Map<String, Integer> oncologyGroupReadCounts,
            double meanReadLength)
    {
        if(meanReadLength <= 0)
        {
            throw new IllegalArgumentException("Invalid mean read length: " + meanReadLength);
        }

        return contigStats.stream()
                .collect(groupingBy(stats -> stats.contig().oncologyGroup()))
                .entrySet().stream()
                .map(entry -> selectOncologyGroup(
                        entry.getKey(), entry.getValue(), margins, readCount(oncologyGroupReadCounts, entry.getKey()),
                        meanReadLength))
                .toList();
    }

    private OncologyGroupSelection selectOncologyGroup(
            String oncologyGroup, List<ContigStats> groupContigs, PairwiseMargins margins, int oncologyGroupReads,
            double meanReadLength)
    {
        GroupCandidates prefiltered = GroupCandidates.prefilter(groupContigs, meanReadLength);
        List<ContigStats> candidates = prefiltered.candidates();

        List<ContigSelectionResult> rejected = prefiltered.rejected().stream()
                .map(rejection -> new ContigSelectionResult(rejection.stats(), rejection.reason(), null))
                .toList();

        if(candidates.isEmpty())
        {
            return new OncologyGroupSelection(oncologyGroup, OncologyGroupOutcome.NO_CANDIDATES, rejected);
        }

        ChallengeGraph graph = ChallengeGraph.build(candidates, oncologyGroupReads, margins);
        RepresentativeChoice choice = RepresentativeChoice.from(candidates, graph);
        Map<ViralContig, Integer> votesRankByContig = votesRanks(candidates);

        List<ContigSelectionResult> candidateResults = candidates.stream()
                .map(candidate -> candidateResult(candidate, choice, graph, candidates, votesRankByContig))
                .toList();

        return new OncologyGroupSelection(
                oncologyGroup, choice.outcome(), Stream.concat(rejected.stream(), candidateResults.stream()).toList());
    }

    private static ContigSelectionResult candidateResult(
            ContigStats candidate, RepresentativeChoice choice, ChallengeGraph graph, List<ContigStats> candidates,
            Map<ViralContig, Integer> votesRankByContig)
    {
        List<Integer> challenges = candidateRanks(
                candidates, candidate, other -> graph.challenges(candidate.contig(), other), votesRankByContig);
        List<Integer> challengedBy = candidateRanks(
                candidates, candidate, other -> graph.challenges(other, candidate.contig()), votesRankByContig);

        CandidateSelectionResult result = new CandidateSelectionResult(
                votesRankByContig.get(candidate.contig()), graph.comparable().contains(candidate.contig()),
                challenges, challengedBy, choice.role(candidate.contig()));

        return new ContigSelectionResult(candidate, ContigFilterStatus.CANDIDATE, result);
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

    private static int readCount(Map<String, Integer> oncologyGroupReadCounts, String oncologyGroup)
    {
        return requireNonNull(
                oncologyGroupReadCounts.get(oncologyGroup), "No aligned read count for oncology group: " + oncologyGroup);
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
}
