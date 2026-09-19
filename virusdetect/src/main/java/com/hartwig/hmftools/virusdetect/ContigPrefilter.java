package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_COVERAGE_LOWER;
import static com.hartwig.hmftools.virusdetect.VirusConstants.MIN_VOTES_PER_BASE;

// Determines whether a contig carries enough evidence to be considered present in the sample.
public class ContigPrefilter
{
    // TODO: this method is kind of weird. Bit off to let the calling code coordinate when this result is just passed back into status(). Could be packaged more cohesively?
    // Any contig above the min coverage establishes that its oncology group is present in the sample.
    public static boolean establishesGroupPresence(double coverageFraction)
    {
        return coverageFraction >= MIN_COVERAGE;
    }

    public static ContigFilterStatus status(
            ViralContig contig, double coverageFraction, double readVotes, boolean groupPresent, double meanReadLength)
    {
        // Once some contig has established the group, its siblings are kept down to a slightly lower coverage, so a
        // near-identical sibling straddling the cutoff is not harshly lost.
        if(!groupPresent || coverageFraction < MIN_COVERAGE_LOWER)
        {
            return ContigFilterStatus.LOW_COVERAGE;
        }
        else if(readVotes < voteFloor(contig, meanReadLength))
        {
            return ContigFilterStatus.LOW_VOTE_DENSITY;
        }
        else
        {
            return ContigFilterStatus.CANDIDATE;
        }
    }

    // Enough votes to have covered the minimum coverage fraction of the contig at the required depth.
    // Drops contigs which are similar to other contigs and so get coverage, but the reads decisively prefer the other
    // contig.
    private static double voteFloor(ViralContig contig, double meanReadLength)
    {
        return MIN_VOTES_PER_BASE * MIN_COVERAGE * contig.length() / meanReadLength;
    }
}
