package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;

import org.jetbrains.annotations.NotNull;

public class ExtractHomozygousDeletes {

    private static final int MIN_WINDOW_COUNT = 10;
    private static final double MAX_COPY_NUMBER = 0.5;
    private static final double MIN_ALLELE_PLOIDY = 0.8;

    @NotNull
    public static List<CombinedRegion> extractHomozygousDeletes(@NotNull final List<CombinedRegion> regions) {
        final List<CombinedRegion> result = Lists.newArrayList();
        for (final CombinedRegion parent : regions) {
            final List<CombinedRegion> children = Lists.newArrayList();
            for (final FittedRegion child : parent.regions()) {
                if (isHomozgousDelete(parent, child)) {
                    children.add(createChild(child));
                }
            }
            result.addAll(ExtractChildren.fillGaps(parent, children));
        }

        return result;
    }

    private static boolean isHomozgousDelete(@NotNull final CombinedRegion parent, @NotNull final FittedRegion child) {
        return child.status() == GermlineStatus.DIPLOID
                && Doubles.lessOrEqual(child.tumorCopyNumber(), MAX_COPY_NUMBER)
                && Doubles.greaterThan(parent.tumorCopyNumber(), MAX_COPY_NUMBER)
                && child.depthWindowCount() >= MIN_WINDOW_COUNT
                && Doubles.greaterOrEqual(parent.region().minorAllelePloidy(), MIN_ALLELE_PLOIDY)
                && Doubles.greaterOrEqual(parent.region().majorAllelePloidy(), MIN_ALLELE_PLOIDY);
    }

    @NotNull
    private static CombinedRegion createChild(@NotNull final FittedRegion child) {
        final CombinedRegion result = new CombinedRegionImpl(child);
        result.setTumorCopyNumber(CopyNumberMethod.HOMOZYGOUS_DELETION, child.tumorCopyNumber());
        return result;
    }

}
