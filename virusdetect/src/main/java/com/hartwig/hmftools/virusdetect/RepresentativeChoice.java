package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.toSet;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.Nullable;

// Which of an oncology group's candidate contigs represents the virus in the sample, decided from the challenges graph,
// with a role for every candidate explaining the decision and an outcome for the group.
// Errs cautious: only clear cases resolve, otherwise the group is left unresolved to avoid silently guessing wrong.
public record RepresentativeChoice(
        Map<ViralContig, ContigRole> roles,
        OncologyGroupOutcome outcome
)
{
    public static RepresentativeChoice from(List<ContigStats> candidates, ChallengeGraph graph)
    {
        List<ViralContig> contigs = graph.contigs();
        if(contigs.size() == 1)
        {
            return new RepresentativeChoice(
                    Map.of(contigs.get(0), ContigRole.REPRESENTATIVE), OncologyGroupOutcome.ONE_CANDIDATE);
        }

        Set<ViralContig> comparable = graph.comparable();

        boolean minorChallengesAbundant = contigs.stream()
                .filter(contig -> !comparable.contains(contig))
                .anyMatch(minor -> comparable.stream().anyMatch(peer -> graph.challenges(minor, peer)));

        Set<ViralContig> challengedByPeer = comparable.stream()
                .filter(contig -> comparable.stream()
                        .anyMatch(peer -> !peer.equals(contig) && graph.challenges(peer, contig)))
                .collect(toSet());
        List<ViralContig> unchallenged = comparable.stream()
                .filter(contig -> !challengedByPeer.contains(contig)).toList();

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
            representative = bestSupported(candidates, unchallenged);
        }

        Map<ViralContig, ContigRole> roles = new HashMap<>();
        for(ViralContig contig : contigs)
        {
            roles.put(contig, roleFor(contig, representative, comparable, challengedByPeer, graph));
        }
        return new RepresentativeChoice(roles, outcome);
    }

    public ContigRole role(ViralContig contig)
    {
        ContigRole role = roles.get(contig);
        if(role == null)
        {
            throw new IllegalArgumentException("contig was not a selection candidate: " + contig.name());
        }
        return role;
    }

    private static ViralContig bestSupported(List<ContigStats> candidates, List<ViralContig> contigs)
    {
        return candidates.stream()
                .filter(stats -> contigs.contains(stats.contig()))
                .min(ContigStats.BEST_SUPPORT_FIRST)
                .orElseThrow()
                .contig();
    }

    private static ContigRole roleFor(
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
        // Abundant and unchallenged by peers, but not crowned: a resolved group makes it a twin, else contested.
        return representative != null ? ContigRole.REPRESENTATIVE_TWIN : ContigRole.CONTESTED;
    }

    private static boolean hasMutualChallenge(Set<ViralContig> comparable, ChallengeGraph graph)
    {
        return comparable.stream().anyMatch(subject -> comparable.stream().anyMatch(opponent ->
                !subject.equals(opponent)
                        && graph.challenges(subject, opponent) && graph.challenges(opponent, subject)));
    }
}
