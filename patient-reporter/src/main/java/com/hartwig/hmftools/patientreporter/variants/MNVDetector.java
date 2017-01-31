package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

final class MNVDetector {

    private static final long MAX_DISTANCE_FOR_MNV = 3;

    private MNVDetector() {
    }

    @NotNull
    static List<SomaticVariant> locatePotentialMNVs(@NotNull final List<SomaticVariant> allVariants,
            @NotNull final List<SomaticVariant> reportedVariants) {
        final List<SomaticVariant> potentialMNVs = Lists.newArrayList();
        for (final SomaticVariant variant : allVariants) {
            reportedVariants.stream().filter(
                    reportVariant -> variant.chromosome().equals(reportVariant.chromosome())).forEach(
                    reportVariant -> {
                        final long distance = Math.abs(variant.position() - reportVariant.position());
                        if (distance > 0 && distance < MAX_DISTANCE_FOR_MNV) {
                            potentialMNVs.add(reportVariant);
                        }
                    });

            // KODU: Performance optimization, once an MNV -> always an MNV!
            potentialMNVs.stream().filter(reportedVariants::contains).forEach(reportedVariants::remove);
        }

        return potentialMNVs;
    }
}
