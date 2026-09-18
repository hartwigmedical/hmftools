package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.COMPARABLE_VOTE_RATIO;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_READS;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

// The directed "challenges" relation over an oncology group's candidate contigs.
// A contig "challenges" another when a high enough fraction of the group's reads fit it better with a high enough margin.
// Also holds which contigs are abundant (vote share near the group's top), the contenders that contest each other.
public class ChallengeGraph
{
    private final List<ViralContig> mContigs;
    private final Set<ViralContig> mComparable;
    // subject contig -> opponent contigs it challenges
    private final Map<ViralContig, Set<ViralContig>> mChallenges;

    private ChallengeGraph(List<ViralContig> contigs, Set<ViralContig> comparable, Map<ViralContig, Set<ViralContig>> challenges)
    {
        mContigs = contigs;
        mComparable = comparable;
        mChallenges = challenges;
    }

    public static ChallengeGraph build(
            List<ContigStats> candidates, Collection<ContigStats> oncologyGroupContigs, PairwiseMargins margins)
    {
        List<ViralContig> contigs = candidates.stream().map(ContigStats::contig).toList();

        double topVotes = candidates.stream().mapToDouble(ContigStats::readVotes).max().orElse(0.0);
        Set<ViralContig> comparable = candidates.stream()
                .filter(stats -> topVotes > 0 && stats.readVotes() >= COMPARABLE_VOTE_RATIO * topVotes)
                .map(ContigStats::contig)
                .collect(toSet());

        double groupVotes = groupVotes(oncologyGroupContigs);
        Map<ViralContig, Set<ViralContig>> challenges = new HashMap<>();
        for(ViralContig subject : contigs)
        {
            challenges.put(
                    subject, contigs.stream()
                            .filter(opponent -> !opponent.equals(subject))
                            .filter(opponent -> challenges(subject, opponent, margins, groupVotes))
                            .collect(toSet()));
        }

        return new ChallengeGraph(contigs, comparable, challenges);
    }

    public List<ViralContig> contigs()
    {
        return mContigs;
    }

    public Set<ViralContig> comparable()
    {
        return mComparable;
    }

    public boolean challenges(ViralContig subject, ViralContig opponent)
    {
        return mChallenges.getOrDefault(subject, Set.of()).contains(opponent);
    }

    private static boolean challenges(ViralContig subject, ViralContig opponent, PairwiseMargins margins, double groupVotes)
    {
        return margins.readsWinningBy(subject, opponent, MIN_CHALLENGE_MARGIN) / groupVotes >= MIN_CHALLENGE_READS;
    }

    private static double groupVotes(Collection<ContigStats> oncologyGroupContigs)
    {
        double groupVotes = oncologyGroupContigs.stream().mapToDouble(ContigStats::readVotes).sum();
        if(groupVotes <= 0)
        {
            String oncologyGroup = oncologyGroupContigs.iterator().next().contig().oncologyGroup();
            throw new IllegalStateException("Oncology group has contig stats but no read votes: " + oncologyGroup);
        }
        return groupVotes;
    }
}
