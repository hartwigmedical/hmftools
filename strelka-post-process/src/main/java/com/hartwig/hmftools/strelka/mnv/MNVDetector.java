package com.hartwig.hmftools.strelka.mnv;

import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class MNVDetector {
    private static final int MNV_GAP_SIZE = 1;

    private MNVDetector() {
    }

    // MIVO: check if variant fits current mnv region. output new region and previous region.
    //  if true, add mnv to current region => return pair of: updated region, empty
    //  if false, we moved past the current region => return pair of: new mnv region starting at variant, previous region (region passed as param)
    @NotNull
    public static Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> addMnvToRegion(
            @NotNull final PotentialMNVRegion potentialMnvRegion, @NotNull final VariantContext variant) {
        return addMnvToRegion(potentialMnvRegion, variant, MNV_GAP_SIZE);
    }

    @VisibleForTesting
    static Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> addMnvToRegion(@NotNull final PotentialMNVRegion potentialMnvRegion,
            @NotNull final VariantContext variant, final int gapSize) {
        if (variantFitsRegion(potentialMnvRegion, variant, gapSize)) {
            return ImmutablePair.of(PotentialMNVRegion.addVariant(potentialMnvRegion, variant, gapSize), Optional.empty());
        } else {
            return ImmutablePair.of(PotentialMNVRegion.fromVariant(variant), Optional.of(potentialMnvRegion));
        }
    }

    @VisibleForTesting
    static boolean variantFitsRegion(@NotNull final PotentialMNVRegion potentialMnvRegion, @NotNull final VariantContext variant,
            final int gapSize) {
        return potentialMnvRegion.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMnvRegion.end() <= gapSize
                && variant.getStart() - potentialMnvRegion.start() >= 0;
    }
}
