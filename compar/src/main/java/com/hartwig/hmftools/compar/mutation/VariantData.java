package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.codon.AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER;
import static com.hartwig.hmftools.common.codon.HgvsCommon.HGVS_STOP_TRI_CODE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_PERC_THRESHOLD;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.SimpleVariant;
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
    public final double AlleleFrequency;
    public final int TumorSupportingReadCount;
    public final int TumorTotalReadCount;

    public final boolean IsFromUnfilteredVcf; // may not use this class for unfiltered variants

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
    static final String FLD_AF = "AF";
    static final String FLD_TUMOR_SUPPORTING_READ_COUNT = "TumorSupportingReadCount";
    static final String FLD_TUMOR_TOTAL_READ_COUNT = "TumorTotalReadCount";

    protected static final double NO_QUAL_PRESENT = -10;

    public VariantData(
            final CategoryType category,
            final String chromosome, final int position, final String ref, final String alt, final VariantType type, final String gene,
            final boolean reported, final HotspotType hotspotStatus, final VariantTier tier,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final int qual,
            final Set<String> filters, final double variantCopyNumber, final double alleleFrequency,
            final int tumorSupportingReadCount, final int tumorTotalReadCount, final boolean isFromUnfilteredVcf)
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
        AlleleFrequency = alleleFrequency;
        TumorSupportingReadCount = tumorSupportingReadCount;
        TumorTotalReadCount = tumorTotalReadCount;

        IsFromUnfilteredVcf = isFromUnfilteredVcf;
    }

    protected void addDefaultValues(final List<FieldInfo> fields)
    {
        addBoolValue(FLD_REPORTED, Reported, fields);
        addStringValue(FLD_HOTSPOT, HotspotStatus.toString(), fields);
        addStringValue(FLD_TIER, Tier.toString(), fields);
        addIntValue(FLD_QUAL, Qual, fields);
        addDoubleValue(FLD_AF, AlleleFrequency, fields);
        addIntValue(FLD_TUMOR_SUPPORTING_READ_COUNT, TumorSupportingReadCount, fields);
        addIntValue(FLD_TUMOR_TOTAL_READ_COUNT, TumorTotalReadCount, fields);
        addStringValue(FLD_FILTER, filtersStr(), fields);

        if(!IsFromUnfilteredVcf)
        {
            // leave out of comparison if created from a filtered variant
            addDoubleValue(FLD_VARIANT_COPY_NUMBER, VariantCopyNumber, fields);
            addStringValue(FLD_GENE, Gene, fields);
            addStringValue(FLD_CANON_EFFECT, CanonicalEffect, fields);
            addStringValue(FLD_CODING_EFFECT, CanonicalCodingEffect, fields);
            addStringValue(FLD_HGVS_CODING, CanonicalHgvsCodingImpact, fields);
            addStringValue(FLD_HGVS_PROTEIN, CanonicalHgvsProteinImpact, fields);
            addStringValue(FLD_OTHER_REPORTED, OtherReportedEffects, fields);
        }
    }

    @Override
    public CategoryType category() { return Category; }

    @Override
    public String key()
    {
        return String.format("%s:%d %s>%s %s", Chromosome, Position, Ref, Alt, Type);
    }

    @Override
    public boolean reportable()
    {
        return !IsFromUnfilteredVcf && Reported;
    }

    @Override
    public boolean isPass()
    {
        return !IsFromUnfilteredVcf && filtersStr().equals(PASS_FILTER);
    }

    @Override
    public String geneName() { return Gene; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final VariantData otherVar = (VariantData) other;

        if(!Chromosome.equals(otherVar.Chromosome) || Position != otherVar.Position)
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
        else
        {
            return PASS_FILTER;
        }
    }

    static List<String> sharedFieldNames()
    {
        return Lists.newArrayList(
                FLD_REPORTED, FLD_HOTSPOT, FLD_TIER, FLD_GENE, FLD_CANON_EFFECT, FLD_CODING_EFFECT,
                FLD_HGVS_CODING, FLD_HGVS_PROTEIN, FLD_OTHER_REPORTED, FLD_QUAL, FLD_VARIANT_COPY_NUMBER, FLD_AF,
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
                getOrMakeFieldCheck(fieldCheckMap, FLD_VARIANT_COPY_NUMBER, COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        fields.add(new FieldInfo(
                FLD_AF,
                getOrMakeFieldCheck(fieldCheckMap, FLD_AF, 0.2, null),
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

    protected static SimpleVariant fromTruthsetKey(final String key)
    {
        try
        {
            String[] keyParts = key.split(":", 4);
            return new SimpleVariant(keyParts[0], Integer.parseInt(keyParts[1]), keyParts[2], keyParts[3]);
        }
        catch(Exception e)
        {
            CMP_LOGGER.error("invalid variant key({})", key);
            return null;
        }
    }

    protected static String checkConvertAminoAcids(final String hgvsProtein)
    {
        // a convenience for truthset variants using single-letter amino acids
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < hgvsProtein.length(); ++i)
        {
            char c = hgvsProtein.charAt(i);
            String cStr = String.valueOf(c);

            // check if the character matches a known single-letter AA
            String longAA = TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.entrySet().stream()
                    .filter(x -> x.getValue().equals(cStr)).map(x -> x.getKey()).findFirst().orElse(null);

            if(longAA == null)
            {
                sb.append(c);
                continue;
            }

            // check if the 3-letter sequence matches a known 3-letter AA
            if(i <= hgvsProtein.length() - 3)
            {
                String existingLongAA = hgvsProtein.substring(i, i + 3);

                if(existingLongAA.equals(longAA)
                || TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.containsKey(existingLongAA)
                || existingLongAA.equals(HGVS_STOP_TRI_CODE))
                {
                    sb.append(existingLongAA);
                    i += 2;
                    continue;
                }
            }

            // if not, use the 3-letter AA`
            sb.append(longAA);
        }

        return sb.toString();
    }
}
