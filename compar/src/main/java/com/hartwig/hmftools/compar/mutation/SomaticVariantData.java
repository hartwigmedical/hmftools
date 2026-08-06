package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_PROB;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.compar.common.CategoryType.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantData extends VariantData
{
    public final boolean Biallelic;
    public final double BiallelicProbability;
    public final boolean HasLPS;
    public final double SubclonalLikelihood;

    protected static final String FLD_BIALLELIC = "Biallelic";
    protected static final String FLD_BIALLELIC_PROB = "BiallelicProb";
    protected static final String FLD_LPS = "HasLPS";
    protected static final String FLD_SUBCLONAL_LIKELIHOOD = "SubclonalLikelihood";

    public SomaticVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type, final String gene,
            final boolean reported, final HotspotType hotspotStatus, final VariantTier tier,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final int qual,
            final Set<String> filters, final double variantCopyNumber, final double purityAdjustedVaf,
            final int tumorSupportingReadCount, final int tumorTotalReadCount, final boolean isFromUnfilteredVcf,
            boolean biallelic, double biallelicProb, final boolean hasLPS, final double subclonalLikelihood,
            final List<FieldInfo> fields)
    {
        super(
                SOMATIC_VARIANT, chromosome, position, ref, alt, type, gene, reported, hotspotStatus, tier,
                canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact, canonicalHgvsProteinImpact, otherReportedEffects,
                qual, filters, variantCopyNumber, purityAdjustedVaf, tumorSupportingReadCount, tumorTotalReadCount,
                isFromUnfilteredVcf, fields);

        BiallelicProbability = biallelicProb;
        Biallelic = biallelic;
        HasLPS = hasLPS;
        SubclonalLikelihood = subclonalLikelihood;
    }

    public static SomaticVariantData fromContext(
            final VariantContext context, final String sampleId, final boolean fromUnfilteredFile,
            final SourceType sourceType, final ComparConfig config, final List<FieldInfo> fields)
    {
        int position = context.getStart();
        String chromosome = context.getContig();
        String ref = context.getReference().getBaseString();
        String alt = !context.getAlternateAlleles().isEmpty() ? context.getAlternateAlleles().get(0).toString() : ref;

        VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(context);

        if(config.RequiresLiftover && sourceType == SourceType.OLD)
        {
            BasePosition comparisonPosition = determineComparisonGenomePosition(
                    chromosome, position, sourceType, config.RequiresLiftover, config.LiftoverCache);

            position = comparisonPosition.Position;
            chromosome = comparisonPosition.Chromosome;
        }

        var tumorAllelicDepth = AllelicDepth.fromGenotype(context.getGenotype(sampleId));

        SomaticVariantData variant = new SomaticVariantData(
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
                context.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false),
                context.getAttributeAsDouble(PURPLE_BIALLELIC_PROB, 0),
                context.hasAttribute(LOCAL_PHASE_SET),
                context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0),
                fields);

        variant.addVAllValues(fields);
        return variant;
    }

    protected void addVAllValues(final List<FieldInfo> fields)
    {
        addDefaultValues(fields);

        addBoolValue(FLD_LPS, HasLPS, fields);

        if(!IsFromUnfilteredVcf)
        {
            addDoubleValue(FLD_BIALLELIC_PROB, BiallelicProbability, fields);
            addBoolValue(FLD_BIALLELIC, Biallelic, fields);
            addDoubleValue(FLD_SUBCLONAL_LIKELIHOOD, SubclonalLikelihood, fields);
        }
    }

    public static SomaticVariantData fromTruthset(final List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        SimpleVariant simpleVariant = fromTruthsetKey(truthsetValues.get(0).Key);

        if(simpleVariant == null)
            return null;

        String chromosome = simpleVariant.Chromosome;
        int position = simpleVariant.Position;
        String ref = simpleVariant.Ref;
        String alt = simpleVariant.Alt;
        VariantType type = simpleVariant.Type;

        String gene = "";
        boolean reported = false;
        HotspotType hotspotStatus = HotspotType.NON_HOTSPOT;
        VariantTier tier = VariantTier.LOW_CONFIDENCE;
        String canonicalEffect = "";
        String canonicalCodingEffect = "";
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        String otherReportedEffects = "";
        int qual = 0;
        Set<String> filters = Sets.newHashSet();
        double variantCopyNumber = 0;
        double purityAdjustedVaf = 0;
        int tumorSupportingReadCount = 0;
        int tumorTotalReadCount = 0;
        boolean biallelic = false;
        double biallelicProb = 0;
        boolean hasLPS = false;
        double subclonalLikelihood= 0;

        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(FLD_REPORTED))
                reported = Boolean.parseBoolean(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(FLD_GENE))
                gene = truthsetValue.Value;
            else if(truthsetValue.FieldName.equals(FLD_PURITY_ADJUSTED_VAF))
                purityAdjustedVaf = Double.parseDouble(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(FLD_VARIANT_COPY_NUMBER))
                variantCopyNumber = Double.parseDouble(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(FLD_CANON_EFFECT))
                canonicalEffect = truthsetValue.Value;
            else if(truthsetValue.FieldName.equals(FLD_TUMOR_TOTAL_READ_COUNT))
                tumorTotalReadCount = Integer.parseInt(truthsetValue.Value);
            else if(truthsetValue.FieldName.equals(FLD_TUMOR_SUPPORTING_READ_COUNT))
                tumorSupportingReadCount = Integer.parseInt(truthsetValue.Value);
        }

        SomaticVariantData variant = new SomaticVariantData(
                chromosome, position, ref, alt, type, gene, reported, hotspotStatus, tier,
                canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact, canonicalHgvsProteinImpact, otherReportedEffects,
                qual, filters, variantCopyNumber, purityAdjustedVaf, tumorSupportingReadCount, tumorTotalReadCount,
                false, biallelic, biallelicProb, hasLPS, subclonalLikelihood, fields);

        variant.addTruthsetValues(truthsetValues, fields);

        return variant;
    }

    private void addTruthsetValues(List<TruthsetValue> truthsetValues, final List<FieldInfo> fields)
    {
        // now add fields
        for(TruthsetValue truthsetValue : truthsetValues)
        {
            if(truthsetValue.FieldName.equals(FLD_REPORTED))
                addBoolValue(FLD_REPORTED, Reported, fields);
            else if(truthsetValue.FieldName.equals(FLD_GENE))
                addStringValue(FLD_GENE, Gene, fields);
            else if(truthsetValue.FieldName.equals(FLD_PURITY_ADJUSTED_VAF))
                addDoubleValue(FLD_PURITY_ADJUSTED_VAF, PurityAdjustedVaf, fields);
            else if(truthsetValue.FieldName.equals(FLD_VARIANT_COPY_NUMBER))
                addDoubleValue(FLD_VARIANT_COPY_NUMBER, VariantCopyNumber, fields);
            else if(truthsetValue.FieldName.equals(FLD_CANON_EFFECT))
                addStringValue(FLD_CANON_EFFECT, CanonicalEffect, fields);
            else if(truthsetValue.FieldName.equals(FLD_TUMOR_TOTAL_READ_COUNT))
                addIntValue(FLD_TUMOR_TOTAL_READ_COUNT, TumorTotalReadCount, fields);
            else if(truthsetValue.FieldName.equals(FLD_TUMOR_SUPPORTING_READ_COUNT))
                addIntValue(FLD_TUMOR_SUPPORTING_READ_COUNT, TumorSupportingReadCount, fields);
        }
    }

    public String toString() { return format("%s gene(%s:%s)", key(), Gene, CanonicalCodingEffect); }
}
