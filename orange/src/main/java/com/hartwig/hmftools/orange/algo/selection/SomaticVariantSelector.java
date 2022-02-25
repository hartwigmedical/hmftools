package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantFactory;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariantSelector {

    private static final Set<String> CUPPA_GENES = Sets.newHashSet("ALB", "SFTPB", "SLC34A2");

    private SomaticVariantSelector() {
    }

    @NotNull
    public static List<ReportableVariant> selectNonDrivers(@NotNull List<SomaticVariant> unreportedVariants,
            @NotNull List<ReportableVariant> reportedSomaticVariants, @NotNull List<ProtectEvidence> evidences) {
        List<ReportableVariant> filtered = Lists.newArrayList();
        for (SomaticVariant variant : unreportedVariants) {
            boolean isHotspot = variant.isHotspot();
            boolean hasEvidence = EvidenceSelector.hasEvidence(evidences, ProtectEventGenerator.variantEvent(variant));
            boolean isCodingAndHasPhasedReportedVariant =
                    !variant.gene().isEmpty() && hasReportedVariantWithPhase(reportedSomaticVariants, variant.topLocalPhaseSet());
            boolean isCuppaRelevantVariant = isRelevantForCuppa(variant);

            if (isHotspot || hasEvidence || isCodingAndHasPhasedReportedVariant || isCuppaRelevantVariant) {
                filtered.add(toReportable(variant));
            }
        }
        return filtered;
    }

    private static boolean hasReportedVariantWithPhase(@NotNull List<ReportableVariant> reportedVariants,
            @Nullable Integer targetPhaseSet) {
        if (targetPhaseSet == null) {
            return false;
        }

        for (ReportableVariant variant : reportedVariants) {
            if (variant.localPhaseSet() != null && variant.localPhaseSet().equals(targetPhaseSet)) {
                return true;
            }
        }

        return false;
    }

    private static boolean isRelevantForCuppa(@NotNull SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && CUPPA_GENES.contains(variant.gene()) && variant.repeatCount() <= 6;
    }

    @NotNull
    private static ReportableVariant toReportable(@NotNull SomaticVariant variant) {
        return ReportableVariantFactory.fromVariant(variant, ReportableVariantSource.SOMATIC).driverLikelihood(Double.NaN).build();
    }
}
