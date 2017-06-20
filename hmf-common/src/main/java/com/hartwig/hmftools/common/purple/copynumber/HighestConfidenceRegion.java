package com.hartwig.hmftools.common.purple.copynumber;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

class HighestConfidenceRegion {

    private final String chromosome;
    private final double purity;

    HighestConfidenceRegion(final String chromosome, final double purity) {
        this.chromosome = chromosome;
        this.purity = purity;
    }

    Optional<PurpleCopyNumber> largestBaf(@NotNull final List<FittedRegion> fittedRegions) {
        return fittedRegions.stream().max(Comparator.comparingInt(ObservedRegion::bafCount)).map(this::best);
    }

    @NotNull
    private PurpleCopyNumber best(@NotNull FittedRegion region) {
        return new PurpleCopyNumberBuilder(purity, region).build();
    }

}
