package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.Doubles;

public class GermlineStatusFactory
{
    private static final double GERMLINE_HOM_DELETION_THRESHOLD = 0.1;
    private static final double GERMLINE_HET_DELETION_THRESHOLD = 0.85;
    private static final double GERMLINE_AMPLIFICATION_THRESHOLD = 1.15;
    private static final double GERMLINE_NOISE_THRESHOLD = 2.2;

    private final CobaltChromosomes mCobaltChromosomes;

    public GermlineStatusFactory(final CobaltChromosomes cobaltChromosomes)
    {
        mCobaltChromosomes = cobaltChromosomes;
    }

    public GermlineStatus calcStatus(final String contig, final double normalRatio, final double tumorRatio, int depthWindowCount)
    {
        if(Doubles.isZero(normalRatio) || !mCobaltChromosomes.contains(contig) || depthWindowCount == 0)
        {
            return UNKNOWN;
        }

        final CobaltChromosome chromosome = mCobaltChromosomes.get(contig);
        double adjustment = chromosome.actualRatio();

        double adjustedHomDeletionThreshold = GERMLINE_HOM_DELETION_THRESHOLD * adjustment;
        if(Doubles.lessThan(normalRatio, adjustedHomDeletionThreshold) && Doubles.lessThan(tumorRatio, adjustedHomDeletionThreshold))
        {
            return HOM_DELETION;
        }

        if(Doubles.lessThan(normalRatio, GERMLINE_HET_DELETION_THRESHOLD * adjustment))
        {
            return HET_DELETION;
        }

        if(Doubles.greaterThan(normalRatio, GERMLINE_NOISE_THRESHOLD * adjustment))
        {
            return NOISE;
        }

        if(Doubles.greaterThan(normalRatio, GERMLINE_AMPLIFICATION_THRESHOLD * adjustment))
        {
            return AMPLIFICATION;
        }

        return DIPLOID;
    }
}
