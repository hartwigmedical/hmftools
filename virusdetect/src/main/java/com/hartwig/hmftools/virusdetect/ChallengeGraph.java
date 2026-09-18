package com.hartwig.hmftools.virusdetect;

import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.virusdetect.VirusConstants.COMPARABLE_VOTE_RATIO;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_MARGIN;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_CHALLENGE_READS;

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
    private final double mTopVoteShare;
    // subject contig -> opponent contigs it challenges
    private final Map<ViralContig, Set<ViralContig>> mChallenges;

    private ChallengeGraph(
            List<ViralContig> contigs, Set<ViralContig> comparable, double topVoteShare,
            Map<ViralContig, Set<ViralContig>> challenges)
    {
        mContigs = contigs;
        mComparable = comparable;
        mTopVoteShare = topVoteShare;
        mChallenges = challenges;
    }

    public static ChallengeGraph build(List<ContigStats> candidates, PairwiseMargins pairwise, double voteTotal)
    {
        List<ViralContig> contigs = candidates.stream().map(ContigStats::contig).toList();

        double topVoteShare = candidates.stream().mapToDouble(stats -> voteShare(stats.readVotes(), voteTotal)).max().orElse(0.0);
        Set<ViralContig> comparable = candidates.stream()
                .filter(stats -> topVoteShare > 0 && voteShare(stats.readVotes(), voteTotal) >= COMPARABLE_VOTE_RATIO * topVoteShare)
                .map(ContigStats::contig)
                .collect(toSet());

        Map<ViralContig, Set<ViralContig>> challenges = new HashMap<>();
        for(ViralContig subject : contigs)
        {
            challenges.put(
                    subject, contigs.stream()
                            .filter(opponent -> !opponent.equals(subject))
                            .filter(opponent -> challenges(subject, opponent, pairwise, voteTotal))
                            .collect(toSet()));
        }

        return new ChallengeGraph(contigs, comparable, topVoteShare, challenges);
    }

    public List<ViralContig> contigs()
    {
        return mContigs;
    }

    public Set<ViralContig> comparable()
    {
        return mComparable;
    }

    public double topVoteShare()
    {
        return mTopVoteShare;
    }

    public boolean challenges(ViralContig subject, ViralContig opponent)
    {
        return mChallenges.getOrDefault(subject, Set.of()).contains(opponent);
    }

    private static boolean challenges(ViralContig subject, ViralContig opponent, PairwiseMargins pairwise, double voteTotal)
    {
        // TODO(AUS-429): the challenge-fraction denominator is voteTotal summed over coverage-surviving contigs only.
        // Whether it should instead span all reads aligning to the oncology group is an open calibration decision.
        if(voteTotal <= 0)
        {
            return false;
        }
        return pairwise.challengeReads(subject, opponent, MIN_CHALLENGE_MARGIN) / voteTotal >= MIN_CHALLENGE_READS;
    }

    private static double voteShare(double votes, double voteTotal)
    {
        return voteTotal > 0 ? votes / voteTotal : 0.0;
    }
}
