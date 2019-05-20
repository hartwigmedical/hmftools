package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_MINOR_ALLELE_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_PLOIDY_INFO;

import com.hartwig.hmftools.common.purple.region.GermlineStatus;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContextBuilder;

public interface PurityAdjustedSomaticVariantBuilder {

    PurityAdjustedSomaticVariantBuilder ploidy(double ploidy);

    PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(double copyNumber);

    PurityAdjustedSomaticVariantBuilder adjustedVAF(double vaf);

    PurityAdjustedSomaticVariantBuilder minorAllelePloidy(double map);

    PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus);

    PurityAdjustedSomaticVariantBuilder biallelic(boolean biallelic);

    @NotNull
    static PurityAdjustedSomaticVariantBuilder fromVariantContextBuilder(@NotNull final VariantContextBuilder builder) {
        return new PurityAdjustedSomaticVariantBuilder() {
            @Override
            public PurityAdjustedSomaticVariantBuilder ploidy(final double ploidy) {
                builder.attribute(PURPLE_PLOIDY_INFO, ploidy);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(final double copyNumber) {
                builder.attribute(PURPLE_CN_INFO, copyNumber);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedVAF(final double vaf) {
                builder.attribute(PURPLE_AF_INFO, vaf);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder minorAllelePloidy(final double map) {
                builder.attribute(PURPLE_MINOR_ALLELE_PLOIDY_INFO, map);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus) {
                builder.attribute(PURPLE_GERMLINE_INFO, germlineStatus.toString());
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder biallelic(final boolean biallelic) {
                builder.attribute(PURPLE_BIALLELIC_FLAG, biallelic);
                return this;
            }
        };
    }
}
