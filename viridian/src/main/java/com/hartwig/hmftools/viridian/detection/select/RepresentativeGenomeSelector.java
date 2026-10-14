package com.hartwig.hmftools.viridian.detection.select;

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

import com.hartwig.hmftools.viridian.detection.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Per oncology group, pick at most one representative virus genome.
// First, virus genomes are filtered on coverage and read votes.
// Then a "challenges" graph determines the presence of "rival" virus genomes - genomes with low overall support,
// but decisively supported by a subset of reads. I.e. one genome in the reference doesn't explain the whole
// virus genome present in the sample.
public class RepresentativeGenomeSelector
{
    public static List<OncologyGroupRepresentativeSelection> select(
            Collection<GenomeSupport> genomeSupport, PairwiseMargins margins, Map<OncologyGroup, Integer> groupReadCounts)
    {
        return genomeSupport.stream()
                .collect(groupingBy(support -> support.genome().oncologyGroup()))
                .entrySet().stream()
                .map(entry -> selectOncologyGroup(entry.getKey(), entry.getValue(), margins, groupReadCounts))
                .toList();
    }

    private static OncologyGroupRepresentativeSelection selectOncologyGroup(
            OncologyGroup oncologyGroup, List<GenomeSupport> groupGenomes, PairwiseMargins margins,
            Map<OncologyGroup, Integer> groupReadCounts)
    {
        List<GenomeSupport> rejected = groupGenomes.stream().filter(support -> !support.isCandidate()).toList();
        List<GenomeSupport> candidates = groupGenomes.stream()
                .filter(GenomeSupport::isCandidate)
                .sorted(comparingDouble(GenomeSupport::readVotes).reversed().thenComparing(GenomeSupport::genome))
                .toList();

        if(candidates.isEmpty())
        {
            return new OncologyGroupRepresentativeSelection(oncologyGroup, OncologyGroupOutcome.NO_CANDIDATES, List.of(), rejected);
        }
        // Note that only groups with enough read coverage reach here.
        // Virus genomes with only origin clipped reads drop out.

        int groupReads = requireNonNull(groupReadCounts.get(oncologyGroup));

        List<VirusGenome> genomes = candidates.stream().map(GenomeSupport::genome).toList();
        Set<VirusGenome> comparable = comparableGenomes(candidates);
        Map<VirusGenome, Set<VirusGenome>> challenges = challenges(genomes, margins, groupReads);
        Map<VirusGenome, Set<VirusGenome>> challengedBy = invertChallengesMap(challenges);

        // A low-abundance virus genome contesting an abundant one is a possible hidden strain, and blocks any verdict.
        Set<VirusGenome> minorChallengers = genomes.stream()
                .filter(genome -> !comparable.contains(genome))
                .filter(genome -> !disjoint(challenges.get(genome), comparable))
                .collect(toSet());

        // Comparable virus genomes no comparable peer challenges. In candidate order, so the first is the best
        // supported.
        List<VirusGenome> leaders = genomes.stream()
                .filter(comparable::contains)
                .filter(genome -> disjoint(challengedBy.get(genome), comparable))
                .toList();

        boolean resolved = minorChallengers.isEmpty() && !leaders.isEmpty();

        List<RepresentativeGenomeCandidate> results = IntStream.range(0, candidates.size())
                .mapToObj(index ->
                {
                    GenomeSupport candidate = candidates.get(index);
                    return new RepresentativeGenomeCandidate(
                            candidate, index + 1, comparable.contains(candidate.genome()), challenges.get(candidate.genome()),
                            challengedBy.get(candidate.genome()),
                            decideGenomeRole(candidate.genome(), comparable, minorChallengers, leaders, resolved));
                })
                .toList();

        return new OncologyGroupRepresentativeSelection(
                oncologyGroup, decideOncologyGroupOutcome(results, comparable, challenges), results, rejected);
    }

    private static GenomeRole decideGenomeRole(
            VirusGenome genome, Set<VirusGenome> comparable, Set<VirusGenome> minorChallengers, List<VirusGenome> leaders,
            boolean resolved)
    {
        if(!comparable.contains(genome))
        {
            return minorChallengers.contains(genome) ? GenomeRole.MINOR_CHALLENGER : GenomeRole.MINOR;
        }
        else if(!leaders.contains(genome))
        {
            return GenomeRole.SECONDARY;
        }
        else if(!resolved)
        {
            return GenomeRole.CONTESTED;
        }
        else
        {
            // The best supported leader is crowned, the rest being indistinguishable from it bar the vote tie-break.
            return genome.equals(leaders.get(0)) ? GenomeRole.REPRESENTATIVE : GenomeRole.REPRESENTATIVE_TWIN;
        }
    }

    private static OncologyGroupOutcome decideOncologyGroupOutcome(
            List<RepresentativeGenomeCandidate> candidates, Set<VirusGenome> comparable, Map<VirusGenome, Set<VirusGenome>> challenges)
    {
        if(anyHasRole(candidates, GenomeRole.MINOR_CHALLENGER))
        {
            return OncologyGroupOutcome.MINOR_RIVAL;
        }
        else if(!anyHasRole(candidates, GenomeRole.REPRESENTATIVE))
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

    private static boolean anyHasRole(List<RepresentativeGenomeCandidate> candidates, GenomeRole role)
    {
        return candidates.stream().anyMatch(candidate -> candidate.role() == role);
    }

    private static Set<VirusGenome> comparableGenomes(List<GenomeSupport> candidates)
    {
        double topVotes = candidates.get(0).readVotes();
        if(topVotes <= 0)
        {
            throw new IllegalStateException("Candidate virus genome has no read votes: " + candidates.get(0).genome());
        }

        return candidates.stream()
                .filter(support -> support.readVotes() >= REPRESENTATIVE_COMPARABLE_VOTE_RATIO * topVotes)
                .map(GenomeSupport::genome)
                .collect(toSet());
    }

    // Subject virus genome -> the opponents it challenges.
    private static Map<VirusGenome, Set<VirusGenome>> challenges(List<VirusGenome> genomes, PairwiseMargins margins, int groupReads)
    {
        return genomes.stream().collect(toMap(
                subject -> subject, subject -> genomes.stream()
                        .filter(opponent -> !opponent.equals(subject))
                        .filter(opponent ->
                                margins.readsWinningBy(subject, opponent, REPRESENTATIVE_CHALLENGE_MARGIN_MIN) / (double) groupReads
                                        >= REPRESENTATIVE_CHALLENGE_READS_MIN)
                        .collect(toSet())));
    }

    private static Map<VirusGenome, Set<VirusGenome>> invertChallengesMap(Map<VirusGenome, Set<VirusGenome>> challenges)
    {
        Map<VirusGenome, Set<VirusGenome>> challengedBy = new HashMap<>();
        challenges.keySet().forEach(genome -> challengedBy.put(genome, new HashSet<>()));
        challenges.forEach((subject, opponents) -> opponents.forEach(opponent -> challengedBy.get(opponent).add(subject)));
        return challengedBy;
    }

    private static boolean hasMutualChallenge(Set<VirusGenome> comparable, Map<VirusGenome, Set<VirusGenome>> challenges)
    {
        return comparable.stream().anyMatch(subject -> challenges.get(subject).stream()
                .filter(comparable::contains)
                .anyMatch(opponent -> challenges.get(opponent).contains(subject)));
    }
}
