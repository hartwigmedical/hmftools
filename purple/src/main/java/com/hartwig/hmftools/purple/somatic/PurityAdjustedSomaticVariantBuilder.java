package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN_INFO;

import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public interface PurityAdjustedSomaticVariantBuilder
{
    PurityAdjustedSomaticVariantBuilder variantCopyNumber(double ploidy);

    PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(double copyNumber);

    PurityAdjustedSomaticVariantBuilder adjustedVAF(double vaf);

    PurityAdjustedSomaticVariantBuilder minorAlleleCopyNumber(double map);

    PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus);

    PurityAdjustedSomaticVariantBuilder biallelic(boolean biallelic);

    @NotNull
    static PurityAdjustedSomaticVariantBuilder fromVariantContext(final VariantContext variantContext)
    {
        return new PurityAdjustedSomaticVariantBuilder()
        {
            @Override
            public PurityAdjustedSomaticVariantBuilder variantCopyNumber(final double ploidy)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_VARIANT_CN_INFO, ploidy);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedCopyNumber(final double copyNumber)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_CN_INFO, copyNumber);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder adjustedVAF(final double vaf)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_AF_INFO, String.format("%.4f", vaf));
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder minorAlleleCopyNumber(final double map)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_MINOR_ALLELE_CN_INFO, map);
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder germlineStatus(@NotNull final GermlineStatus germlineStatus)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_GERMLINE_INFO, germlineStatus.toString());
                return this;
            }

            @Override
            public PurityAdjustedSomaticVariantBuilder biallelic(final boolean biallelic)
            {
                variantContext.getCommonInfo().putAttribute(PURPLE_BIALLELIC_FLAG, biallelic);
                return this;
            }
        };
    }
}
