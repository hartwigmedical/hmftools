package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.combine.CombinedRegion;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

final class PopulateUnknown {

    private final Gender gender;

    PopulateUnknown(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    List<CombinedRegion> populateUnknown(@NotNull final List<CombinedRegion> regions) {

        for (CombinedRegion region : regions) {
            if (region.copyNumberMethod() == CopyNumberMethod.UNKNOWN) {
                region.setTumorCopyNumber(CopyNumberMethod.UNKNOWN,
                        HumanChromosome.fromString(region.chromosome()).isDiploid(gender) ? 2 : 1);
            }
        }
        return regions;
    }

}