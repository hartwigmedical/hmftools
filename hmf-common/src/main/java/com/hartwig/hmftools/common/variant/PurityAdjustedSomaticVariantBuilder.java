package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantHeader.PURPLE_VARIANT_CN_INFO;

import com.hartwig.hmftools.common.purple.region.GermlineStatus;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public interface PurityAdjustedSomaticVariantBuilder {

    PurityAdjustedSomaticVariantBuilder ploidy(double ploidy);

    PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(double copyNumber);

    PurityAdjustedSomaticVariantBuilder adjustedVAF(double vaf);

    PurityAdjustedSomaticVariantBuilder minorAllelePloidy(double map);

    PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus);

    PurityAdjustedSomaticVariantBuilder biallelic(boolean biallelic);

    @NotNull
    static PurityAdjustedSomaticVariantBuilder fromVariantContex(@NotNull final VariantContext builder) {
        return new PurityAdjustedSomaticVariantBuilder() {
            @Override
            public PurityAdjustedSomaticVariantBuilder ploidy(final double ploidy) {
                builder.getCommonInfo().putAttribute(PURPLE_VARIANT_CN_INFO, ploidy);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(final double copyNumber) {
                builder.getCommonInfo().putAttribute(PURPLE_CN_INFO, copyNumber);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedVAF(final double vaf) {
                builder.getCommonInfo().putAttribute(PURPLE_AF_INFO, String.format("%.4f", vaf));
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder minorAllelePloidy(final double map) {
                builder.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, map);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus) {
                builder.getCommonInfo().putAttribute(PURPLE_GERMLINE_INFO, germlineStatus.toString());
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder biallelic(final boolean biallelic) {
                builder.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, biallelic);
                return this;
            }
        };
    }
}
