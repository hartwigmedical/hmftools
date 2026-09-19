package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE_LOWER;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_VOTES_PER_BASE;

import java.util.ArrayList;
import java.util.List;

// An oncology group's contigs split by whether they carry enough evidence to proceed to representative selection.
// Rejected contigs are those filtered up-front because they have too little support to be plausible.
public record GroupCandidates(
        List<ContigSupport> candidates,
        List<RejectedContig> rejected
)
{
    public static GroupCandidates prefilter(List<ContigSupport> oncologyGroupContigs, double meanReadLength)
    {
        boolean groupPresent = oncologyGroupContigs.stream()
                .anyMatch(stats -> stats.coverageFraction() >= MIN_COVERAGE);

        List<ContigSupport> candidates = new ArrayList<>();
        List<RejectedContig> rejected = new ArrayList<>();
        for(ContigSupport stats : oncologyGroupContigs)
        {
            ContigFilterStatus status = filterStatus(stats, groupPresent, meanReadLength);
            if(status == ContigFilterStatus.CANDIDATE)
            {
                candidates.add(stats);
            }
            else
            {
                rejected.add(new RejectedContig(stats, status));
            }
        }

        return new GroupCandidates(candidates, rejected);
    }

    private static ContigFilterStatus filterStatus(ContigSupport stats, boolean groupPresent, double meanReadLength)
    {
        // A group counts as present only once some contig clears the coverage floor. Its remaining contigs are kept
        // down to a slightly lower floor, so a near-identical sibling straddling the cutoff is not lost.
        if(!groupPresent || stats.coverageFraction() < MIN_COVERAGE_LOWER)
        {
            return ContigFilterStatus.LOW_COVERAGE;
        }
        if(!passesVoteDensity(stats, meanReadLength))
        {
            return ContigFilterStatus.LOW_VOTE_DENSITY;
        }
        return ContigFilterStatus.CANDIDATE;
    }

    // Drops a contig with enough coverage but few votes, its reads aligning better to a different contig.
    private static boolean passesVoteDensity(ContigSupport stats, double meanReadLength)
    {
        double voteFloorPerBase = MIN_VOTES_PER_BASE * MIN_COVERAGE / meanReadLength;
        return stats.readVotes() >= voteFloorPerBase * stats.contig().length();
    }
}
