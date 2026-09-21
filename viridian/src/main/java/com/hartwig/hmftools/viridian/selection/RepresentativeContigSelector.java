package com.hartwig.hmftools.viridian.selection;

import static java.util.Collections.disjoint;
import static java.util.Comparator.comparingDouble;
import static java.util.Objects.requireNonNull;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_CHALLENGE_MARGIN_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_CHALLENGE_READS_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.REPRESENTATIVE_COMPARABLE_VOTE_RATIO;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import com.hartwig.hmftools.viridian.detection.contig_support.ContigSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Per oncology group, pick at most one representative contig.
// First, contigs are filtered on coverage and read votes.
// Then a "challenges" graph determines the presence of "rival" contigs - contigs with low overall support, but
// decisively supported by a subset of reads. I.e. 1 contig doesn't explain the whole viral genome in the sample.
public class RepresentativeContigSelector
{
    public static List<OncologyGroupRepresentativeSelection> select(
            Collection<ContigSupport> contigSupport, PairwiseMargins margins, Map<OncologyGroup, Integer> groupReadCounts)
    {
        return contigSupport.stream()
                .collect(groupingBy(support -> support.contig().oncologyGroup()))
                .entrySet().stream()
                .map(entry -> selectOncologyGroup(entry.getKey(), entry.getValue(), margins, groupReadCounts))
                .toList();
    }

    private static OncologyGroupRepresentativeSelection selectOncologyGroup(
            OncologyGroup oncologyGroup, List<ContigSupport> groupContigs, PairwiseMargins margins,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        List<ContigSupport> rejected = groupContigs.stream().filter(support -> !support.isCandidate()).toList();
        List<ContigSupport> candidates = groupContigs.stream()
                .filter(ContigSupport::isCandidate)
                .sorted(comparingDouble(ContigSupport::readVotes).reversed().thenComparing(ContigSupport::contig))
                .toList();

        if(candidates.isEmpty())
        {
            return new OncologyGroupRepresentativeSelection(oncologyGroup, OncologyGroupOutcome.NO_CANDIDATES, List.of(), rejected);
        }
        // Note that only groups with enough read coverage reach here.
        // Contigs with only origin clipped reads drop out.

        int groupReads = requireNonNull(groupReadCounts.get(oncologyGroup));

        List<ViralContig> contigs = candidates.stream().map(ContigSupport::contig).toList();
        Set<ViralContig> comparable = comparableContigs(candidates);
        Map<ViralContig, Set<ViralContig>> challenges = challenges(contigs, margins, groupReads);
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

        boolean resolved = minorChallengers.isEmpty() && !leaders.isEmpty();

        List<RepresentativeContigCandidate> results = IntStream.range(0, candidates.size())
                .mapToObj(index ->
                {
                    ContigSupport candidate = candidates.get(index);
                    return new RepresentativeContigCandidate(
                            candidate, index + 1, comparable.contains(candidate.contig()), challenges.get(candidate.contig()),
                            challengedBy.get(candidate.contig()),
                            decideContigRole(candidate.contig(), comparable, minorChallengers, leaders, resolved));
                })
                .toList();

        return new OncologyGroupRepresentativeSelection(
                oncologyGroup, decideOncologyGroupOutcome(results, comparable, challenges), results, rejected);
    }

    private static ContigRole decideContigRole(
            ViralContig contig, Set<ViralContig> comparable, Set<ViralContig> minorChallengers, List<ViralContig> leaders,
            boolean resolved)
    {
        if(!comparable.contains(contig))
        {
            return minorChallengers.contains(contig) ? ContigRole.MINOR_CHALLENGER : ContigRole.MINOR;
        }
        else if(!leaders.contains(contig))
        {
            return ContigRole.SECONDARY;
        }
        else if(!resolved)
        {
            return ContigRole.CONTESTED;
        }
        else
        {
            // The best supported leader is crowned, the rest being indistinguishable from it bar the vote tie-break.
            return contig.equals(leaders.get(0)) ? ContigRole.REPRESENTATIVE : ContigRole.REPRESENTATIVE_TWIN;
        }
    }

    private static OncologyGroupOutcome decideOncologyGroupOutcome(
            List<RepresentativeContigCandidate> candidates, Set<ViralContig> comparable, Map<ViralContig, Set<ViralContig>> challenges)
    {
        if(anyHasRole(candidates, ContigRole.MINOR_CHALLENGER))
        {
            return OncologyGroupOutcome.MINOR_RIVAL;
        }
        else if(!anyHasRole(candidates, ContigRole.REPRESENTATIVE))
        {
            if(hasMutualChallenge(comparable, challenges))
            {
                return OncologyGroupOutcome.MUTUAL;
            }
            else
            {
                return OncologyGroupOutcome.CYCLE;
            }
        }
        else if(candidates.size() == 1)
        {
            return OncologyGroupOutcome.ONE_CANDIDATE;
        }
        else
        {
            return OncologyGroupOutcome.RESOLVED_CANDIDATES;
        }
    }

    private static boolean anyHasRole(List<RepresentativeContigCandidate> candidates, ContigRole role)
    {
        return candidates.stream().anyMatch(candidate -> candidate.role() == role);
    }

    private static Set<ViralContig> comparableContigs(List<ContigSupport> candidates)
    {
        double topVotes = candidates.get(0).readVotes();
        if(topVotes <= 0)
        {
            throw new IllegalStateException("Candidate contig has no read votes: " + candidates.get(0).contig());
        }

        return candidates.stream()
                .filter(support -> support.readVotes() >= REPRESENTATIVE_COMPARABLE_VOTE_RATIO * topVotes)
                .map(ContigSupport::contig)
                .collect(toSet());
    }

    // Subject contig -> the opponents it challenges.
    private static Map<ViralContig, Set<ViralContig>> challenges(List<ViralContig> contigs, PairwiseMargins margins, int groupReads)
    {
        return contigs.stream().collect(toMap(
                subject -> subject, subject -> contigs.stream()
                        .filter(opponent -> !opponent.equals(subject))
                        .filter(opponent ->
                                margins.readsWinningBy(subject, opponent, REPRESENTATIVE_CHALLENGE_MARGIN_MIN) / (double) groupReads
                                        >= REPRESENTATIVE_CHALLENGE_READS_MIN)
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
