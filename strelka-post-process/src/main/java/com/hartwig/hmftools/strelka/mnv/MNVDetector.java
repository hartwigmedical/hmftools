package com.hartwig.hmftools.strelka.mnv;

import java.util.Optional;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class MNVDetector {
    private MNVDetector() {
    }

    // MIVO: check if variant fits current mnv region. output new region and previous region.
    //  if true, add mnv to current region => return pair of: updated region, empty
    //  if false, we moved past the current region => return pair of: new mnv region starting at variant, previous region (region passed as param)
    @NotNull
    public static Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> fitsMNVRegion(@NotNull final PotentialMNVRegion potentialMnvRegion,
            @NotNull final VariantContext variant) {
        if (potentialMnvRegion.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMnvRegion.end() <= 1
                && variant.getStart() - potentialMnvRegion.start() >= 0) {
            return ImmutablePair.of(PotentialMNVRegion.addVariant(potentialMnvRegion, variant), Optional.empty());
        } else {
            return ImmutablePair.of(PotentialMNVRegion.fromVariant(variant), Optional.of(potentialMnvRegion));
        }
    }
}
