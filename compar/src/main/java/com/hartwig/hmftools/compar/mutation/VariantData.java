package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class VariantData extends ComparableItem
{
    public final CategoryType Category;

    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;
    public final String Gene;

    public final boolean Reported;
    public final HotspotType HotspotStatus;
    public final VariantTier Tier;
    public final String CanonicalEffect;
    public final String CanonicalCodingEffect;
    public final String CanonicalHgvsCodingImpact;
    public final String CanonicalHgvsProteinImpact;
    public final String OtherReportedEffects;
    public final int Qual;
    public final Set<String> Filters;
    public final double VariantCopyNumber;
    public final double PurityAdjustedVaf;
    public final int TumorSupportingReadCount;
    public final int TumorTotalReadCount;

    public final boolean IsFromUnfilteredVcf; // may not use this class for unfiltered variants

    public final String mComparisonChromosome;
    public final int mComparisonPosition;

    static final String FLD_REPORTED = "Reported";
    static final String FLD_QUAL = "Qual";
    static final String FLD_FILTER = "Filter";
    static final String FLD_HOTSPOT = "Hotspot";
    static final String FLD_TIER = "Tier";
    static final String FLD_GENE = "Gene";
    static final String FLD_CANON_EFFECT = "CanonicalEffect";
    static final String FLD_CODING_EFFECT = "CanonicalCodingEffect";
    static final String FLD_HGVS_CODING = "CanonicalHgvsCoding";
    static final String FLD_HGVS_PROTEIN = "CanonicalHgvsProtein";
    static final String FLD_OTHER_REPORTED = "OtherReportedEffects";
    static final String FLD_VARIANT_COPY_NUMBER = "VariantCopyNumber";
    static final String FLD_PURITY_ADJUSTED_VAF = "PurityAdjustedVaf";
    static final String FLD_TUMOR_SUPPORTING_READ_COUNT = "TumorSupportingReadCount";
    static final String FLD_TUMOR_TOTAL_READ_COUNT = "TumorTotalReadCount";

    protected static final double NO_QUAL_PRESENT = -10;

    public VariantData(
            final CategoryType category,
            final String chromosome, final int position, final String ref, final String alt, final VariantType type, final String gene,
            final boolean reported, final HotspotType hotspotStatus, final VariantTier tier,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final int qual,
            final Set<String> filters, final double variantCopyNumber, final double purityAdjustedVaf,
            final int tumorSupportingReadCount, final int tumorTotalReadCount, final boolean isFromUnfilteredVcf,
            final String comparisonChromosome, final int comparisonPosition, final List<FieldInfo> fields)
    {
        Category = category;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;
        Reported = reported;
        HotspotStatus = hotspotStatus;
        Tier = tier;
        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        CanonicalHgvsCodingImpact = canonicalHgvsCodingImpact;
        CanonicalHgvsProteinImpact = canonicalHgvsProteinImpact;
        OtherReportedEffects = otherReportedEffects != null ? otherReportedEffects : "";
        Qual = qual;
        Filters = filters;
        VariantCopyNumber = variantCopyNumber;
        PurityAdjustedVaf = purityAdjustedVaf;
        TumorSupportingReadCount = tumorSupportingReadCount;
        TumorTotalReadCount = tumorTotalReadCount;

        IsFromUnfilteredVcf = isFromUnfilteredVcf;

        mComparisonChromosome = comparisonChromosome;
        mComparisonPosition = comparisonPosition;

        addBoolValue(FLD_REPORTED, reported, fields);
        addStringValue(FLD_HOTSPOT, hotspotStatus.toString(), fields);
        addStringValue(FLD_TIER, tier.toString(), fields);
        addStringValue(FLD_GENE, gene, fields);
        addStringValue(FLD_CANON_EFFECT, canonicalEffect, fields);
        addStringValue(FLD_CODING_EFFECT, canonicalCodingEffect, fields);
        addStringValue(FLD_HGVS_CODING, canonicalHgvsCodingImpact, fields);
        addStringValue(FLD_HGVS_PROTEIN, canonicalHgvsProteinImpact, fields);
        addStringValue(FLD_OTHER_REPORTED, otherReportedEffects, fields);

        addIntValue(FLD_QUAL, qual, fields);

        addDoubleValue(FLD_VARIANT_COPY_NUMBER, variantCopyNumber, fields);
        addDoubleValue(FLD_PURITY_ADJUSTED_VAF, purityAdjustedVaf, fields);

        addIntValue(FLD_TUMOR_SUPPORTING_READ_COUNT, tumorSupportingReadCount, fields);
        addIntValue(FLD_TUMOR_TOTAL_READ_COUNT, tumorTotalReadCount, fields);

        addStringValue(FLD_FILTER, filtersStr(), fields);
    }

    @Override
    public CategoryType category() { return Category; }

    @Override
    public String key()
    {
        if(mComparisonPosition != Position)
        {
            return String.format("%s:%d %s>%s %s liftover(%s)",
                    Chromosome, Position, Ref, Alt, Type, mComparisonPosition);
        }
        else
        {
            return String.format("%s:%d %s>%s %s", Chromosome, Position, Ref, Alt, Type);
        }
    }

    @Override
    public boolean reportable()
    {
        return !IsFromUnfilteredVcf && Reported;
    }

    @Override
    public boolean isPass()
    {
        // a reportable variant not in a gene should be impossible
        return !IsFromUnfilteredVcf && (Reported || !Gene.isEmpty());
    }

    @Override
    public String geneName() { return Gene; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final VariantData otherVar = (VariantData) other;

        if(!mComparisonChromosome.equals(otherVar.Chromosome) || mComparisonPosition != otherVar.Position)
            return false;

        if(!Ref.equals(otherVar.Ref) || !Alt.equals(otherVar.Alt))
            return false;

        if(Type != otherVar.Type)
            return false;

        return true;
    }

    public String filtersStr()
    {
        if(!Filters.isEmpty())
        {
            return Filters.stream().collect(Collectors.joining(ITEM_DELIM));
        }
        /*
        else if(IsFromUnfilteredVcf)
        {
            return "FILTERED";
        }
        */
        else
        {
            return PASS_FILTER;
        }
    }

    public String comparisonChromosome() { return mComparisonChromosome; }
    public int comparisonPosition() { return mComparisonPosition; }

    static List<String> sharedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL, FLD_VARIANT_COPY_NUMBER, FLD_PURITY_ADJUSTED_VAF,
                FLD_TUMOR_SUPPORTING_READ_COUNT, FLD_TUMOR_TOTAL_READ_COUNT, FLD_FILTER);
    }

    static void addComparerFields(final List<FieldInfo> fields, final Map<String,FieldCheck> fieldCheckMap)
    {
        fields.add(new FieldInfo(FLD_REPORTED, getOrMakeFieldCheck(fieldCheckMap, FLD_REPORTED), null));
        fields.add(new FieldInfo(FLD_HOTSPOT, getOrMakeFieldCheck(fieldCheckMap, FLD_HOTSPOT), null));
        fields.add(new FieldInfo(FLD_TIER, getOrMakeFieldCheck(fieldCheckMap, FLD_TIER), null));
        fields.add(new FieldInfo(FLD_GENE, getOrMakeFieldCheck(fieldCheckMap, FLD_GENE), null));
        fields.add(new FieldInfo(FLD_CANON_EFFECT, getOrMakeFieldCheck(fieldCheckMap, FLD_CANON_EFFECT), null));
        fields.add(new FieldInfo(FLD_CODING_EFFECT, getOrMakeFieldCheck(fieldCheckMap, FLD_CODING_EFFECT), null));
        fields.add(new FieldInfo(FLD_HGVS_CODING, getOrMakeFieldCheck(fieldCheckMap, FLD_HGVS_CODING), null));
        fields.add(new FieldInfo(FLD_HGVS_PROTEIN, getOrMakeFieldCheck(fieldCheckMap, FLD_HGVS_PROTEIN), null));
        fields.add(new FieldInfo(FLD_OTHER_REPORTED, getOrMakeFieldCheck(fieldCheckMap, FLD_OTHER_REPORTED), null));

        fields.add(new FieldInfo(
                FLD_QUAL,
                getOrMakeFieldCheck(fieldCheckMap, FLD_QUAL, 50.0, 0.2),
                null));

        fields.add(new FieldInfo(
                FLD_VARIANT_COPY_NUMBER,
                getOrMakeFieldCheck(fieldCheckMap, FLD_VARIANT_COPY_NUMBER, 0.3, 0.3),
                "%.2f"));

        fields.add(new FieldInfo(
                FLD_PURITY_ADJUSTED_VAF,
                getOrMakeFieldCheck(fieldCheckMap, FLD_PURITY_ADJUSTED_VAF, 0.2, null),
                "%.2f"));

        fields.add(new FieldInfo(
                FLD_TUMOR_SUPPORTING_READ_COUNT,
                getOrMakeFieldCheck(fieldCheckMap, FLD_TUMOR_SUPPORTING_READ_COUNT, 1.0, 0.2),
                null));

        fields.add(new FieldInfo(
                FLD_TUMOR_TOTAL_READ_COUNT,
                getOrMakeFieldCheck(fieldCheckMap, FLD_TUMOR_TOTAL_READ_COUNT, 1.0, 0.2),
                null));

        fields.add(new FieldInfo(FLD_FILTER, getOrMakeFieldCheck(fieldCheckMap, FLD_FILTER), null));
    }
}
