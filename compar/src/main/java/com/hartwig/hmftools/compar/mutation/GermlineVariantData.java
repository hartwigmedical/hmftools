package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineVariantData extends VariantData
{
    public GermlineVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type, final String gene,
            final boolean reported, final HotspotType hotspotStatus, final VariantTier tier,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final int qual,
            final Set<String> filters, final double variantCopyNumber, final double purityAdjustedVaf,
            final int tumorSupportingReadCount, final int tumorTotalReadCount, final boolean isFromUnfilteredVcf,
            final String comparisonChromosome, final int comparisonPosition,
            final List<FieldInfo> fields)
    {
        super(
                GERMLINE_VARIANT, chromosome, position, ref, alt, type, gene, reported, hotspotStatus, tier,
                canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact, canonicalHgvsProteinImpact, otherReportedEffects,
                qual, filters, variantCopyNumber, purityAdjustedVaf, tumorSupportingReadCount, tumorTotalReadCount,
                isFromUnfilteredVcf, comparisonChromosome, comparisonPosition, fields);
    }

    public static GermlineVariantData fromContext(
            final VariantContext context, final String sampleId, final boolean fromUnfilteredFile,
            final SourceType sourceType, final ComparConfig config, final List<FieldInfo> fields)
    {
        int position = context.getStart();
        String chromosome = context.getContig();
        String ref = context.getReference().getBaseString();
        String alt = !context.getAlternateAlleles().isEmpty() ? context.getAlternateAlleles().get(0).toString() : ref;

        VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(context);

        BasePosition comparisonPosition = determineComparisonGenomePosition(
                chromosome, position, sourceType, config.RequiresLiftover, config.LiftoverCache);

        var tumorAllelicDepth = AllelicDepth.fromGenotype(context.getGenotype(sampleId));

        return new GermlineVariantData(
                chromosome, position, ref, alt, VariantType.type(context),
                variantImpact.GeneName,
                context.getAttributeAsBoolean(REPORTED_FLAG, false),
                HotspotType.fromVariant(context),
                VariantTier.fromContext(context),
                variantImpact.CanonicalEffect,
                variantImpact.CanonicalCodingEffect.toString(),
                variantImpact.CanonicalHgvsCoding,
                variantImpact.CanonicalHgvsProtein,
                variantImpact.OtherReportableEffects,
                (int)context.getPhredScaledQual(),
                context.getFilters(),
                context.getAttributeAsDouble(PURPLE_VARIANT_CN, 0),
                context.getAttributeAsDouble(PURPLE_AF, 0),
                tumorAllelicDepth.AlleleReadCount,
                tumorAllelicDepth.TotalReadCount,
                fromUnfilteredFile,
                comparisonPosition.Chromosome,
                comparisonPosition.Position,
                fields);
    }
}
