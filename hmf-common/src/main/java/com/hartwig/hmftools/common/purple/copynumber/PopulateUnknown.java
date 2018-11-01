package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

final class PopulateUnknown {

    private final Gender gender;

    PopulateUnknown(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    List<CombinedRegion> populateUnknown(@NotNull final List<CombinedRegion> regions) {

        for (int i = 0; i < regions.size(); i++) {
            final CombinedRegion region = regions.get(i);

            if (region.copyNumberMethod() == CopyNumberMethod.UNKNOWN) {

                double normalCopyNumber = HumanChromosome.fromString(region.chromosome()).isDiploid(gender) ? 2 : 1;
                region.setTumorCopyNumber(CopyNumberMethod.UNKNOWN, normalCopyNumber);

                if (region.support() == SegmentSupport.NONE && i > 0) {
                    final CombinedRegion prev = regions.get(i - 1);
                    if (prev.copyNumberMethod() == CopyNumberMethod.UNKNOWN) {
                        prev.extend(region.region());
                        regions.remove(i);
                        i--;
                    }
                }
            }

        }
        return regions;
    }

}