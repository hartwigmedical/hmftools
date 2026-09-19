package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.disjoint;
import static java.util.Comparator.comparingDouble;
import static java.util.Objects.requireNonNull;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.COMPARABLE_VOTE_RATIO;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_READS;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.Nullable;

// Per oncology group, pick at most one representative contig.
// First, contigs are filtered on coverage and read votes.
// Then a "challenges" graph determines the presence of "rival" contigs - contigs with low overall support, but
// decisively supported by a subset of reads. I.e. 1 contig doesn't explain the whole viral genome in the sample.
public class RepresentativeSelector
{
    public List<OncologyGroupSelection> select(
            Collection<ContigSupport> contigSupport, PairwiseMargins margins, Map<OncologyGroup, Integer> groupReadCounts)
    {
        return contigSupport.stream()
                .collect(groupingBy(support -> support.contig().oncologyGroup()))
                .entrySet().stream()
                .map(entry -> selectOncologyGroup(entry.getKey(), entry.getValue(), margins, groupReadCounts))
                .toList();
    }

    private static OncologyGroupSelection selectOncologyGroup(
            OncologyGroup oncologyGroup, List<ContigSupport> groupContigs, PairwiseMargins margins,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        List<ContigSupport> rejected = groupContigs.stream().filter(support -> !support.isCandidate()).toList();
        List<ContigSupport> candidates = groupContigs.stream()
                .filter(ContigSupport::isCandidate)
                .sorted(comparingDouble(ContigSupport::readVotes).reversed().thenComparing(support -> support.contig().name()))
                .toList();

        if(candidates.isEmpty())
        {
            return new OncologyGroupSelection(oncologyGroup, OncologyGroupOutcome.NO_CANDIDATES, List.of(), rejected);
        }

        List<ViralContig> contigs = candidates.stream().map(ContigSupport::contig).toList();
        Set<ViralContig> comparable = comparableContigs(candidates);
        Map<ViralContig, Set<ViralContig>> challenges = challenges(contigs, margins, oncologyGroup, groupReadCounts);
        Map<ViralContig, Set<ViralContig>> challengedBy = invertChallengesMap(challenges);

        // A low-abundance contig contesting an abundant one is a possible hidden strain, and blocks any verdict.
        Set<ViralContig> minorChallengers = contigs.stream()
                .filter(contig -> !comparable.contains(contig))
                .filter(contig -> !disjoint(challenges.get(contig), comparable))
                .collect(toSet());

        // Comparable contigs no comparable peer challenges. In candidate order, so the first is the best supported.
        List<ViralContig> leaders = contigs.stream()
                .filter(comparable::contains)
                .filter(contig -> disjoint(challengedBy.get(contig), comparable))
                .toList();

        OncologyGroupOutcome outcome = decideOncologyGroupOutcome(candidates.size(), minorChallengers, leaders, comparable, challenges);
        ViralContig representative = outcome.resolution() == OncologyGroupResolution.RESOLVED ? leaders.get(0) : null;
        Set<ViralContig> leadSet = Set.copyOf(leaders);

        List<RepresentativeCandidate> results = candidates.stream()
                .map(candidate -> new RepresentativeCandidate(
                        candidate, comparable.contains(candidate.contig()), challenges.get(candidate.contig()),
                        challengedBy.get(candidate.contig()),
                        decideContigRole(candidate.contig(), representative, comparable, minorChallengers, leadSet)))
                .toList();

        return new OncologyGroupSelection(oncologyGroup, outcome, results, rejected);
    }

    private static OncologyGroupOutcome decideOncologyGroupOutcome(
            int candidateCount, Set<ViralContig> minorChallengers, List<ViralContig> leaders, Set<ViralContig> comparable,
            Map<ViralContig, Set<ViralContig>> challenges)
    {
        if(!minorChallengers.isEmpty())
        {
            return OncologyGroupOutcome.MINOR_RIVAL;
        }
        else if(leaders.isEmpty())
        {
            if(hasMutualChallenge(comparable, challenges))
            {
                return OncologyGroupOutcome.MUTUAL;
            }
            return OncologyGroupOutcome.CYCLE;
        }
        else if(candidateCount == 1)
        {
            return OncologyGroupOutcome.ONE_CANDIDATE;
        }
        else
        {
            return OncologyGroupOutcome.RESOLVED_CANDIDATES;
        }
    }

    private static ContigRole decideContigRole(
            ViralContig contig, @Nullable ViralContig representative, Set<ViralContig> comparable,
            Set<ViralContig> minorChallengers, Set<ViralContig> leaders)
    {
        if(contig.equals(representative))
        {
            return ContigRole.REPRESENTATIVE;
        }
        if(!comparable.contains(contig))
        {
            return minorChallengers.contains(contig) ? ContigRole.MINOR_CHALLENGER : ContigRole.MINOR;
        }
        if(!leaders.contains(contig))
        {
            return ContigRole.SECONDARY;
        }
        // Abundant and unchallenged by peers, but not crowned: a resolved group makes it a twin, else contested.
        return representative != null ? ContigRole.REPRESENTATIVE_TWIN : ContigRole.CONTESTED;
    }

    private static Set<ViralContig> comparableContigs(List<ContigSupport> candidates)
    {
        double topVotes = candidates.get(0).readVotes();
        if(topVotes <= 0)
        {
            throw new IllegalStateException("Candidate contig has no read votes: " + candidates.get(0).contig().name());
        }

        return candidates.stream()
                .filter(support -> support.readVotes() >= COMPARABLE_VOTE_RATIO * topVotes)
                .map(ContigSupport::contig)
                .collect(toSet());
    }

    // Subject contig -> the opponents it challenges.
    private static Map<ViralContig, Set<ViralContig>> challenges(
            List<ViralContig> contigs, PairwiseMargins margins, OncologyGroup oncologyGroup,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        int groupReads = requireNonNull(
                groupReadCounts.get(oncologyGroup), "No aligned read count for oncology group: " + oncologyGroup);

        return contigs.stream().collect(toMap(
                subject -> subject, subject -> contigs.stream()
                        .filter(opponent -> !opponent.equals(subject))
                        .filter(opponent -> margins.readsWinningBy(subject, opponent, MIN_CHALLENGE_MARGIN) / (double) groupReads
                                >= MIN_CHALLENGE_READS)
                        .collect(toSet())));
    }

    private static Map<ViralContig, Set<ViralContig>> invertChallengesMap(Map<ViralContig, Set<ViralContig>> challenges)
    {
        Map<ViralContig, Set<ViralContig>> challengedBy = new HashMap<>();
        challenges.keySet().forEach(contig -> challengedBy.put(contig, new HashSet<>()));
        challenges.forEach((subject, opponents) -> opponents.forEach(opponent -> challengedBy.get(opponent).add(subject)));
        return challengedBy;
    }

    private static boolean hasMutualChallenge(Set<ViralContig> comparable, Map<ViralContig, Set<ViralContig>> challenges)
    {
        return comparable.stream().anyMatch(subject -> challenges.get(subject).stream()
                .filter(comparable::contains)
                .anyMatch(opponent -> challenges.get(opponent).contains(subject)));
    }
}
