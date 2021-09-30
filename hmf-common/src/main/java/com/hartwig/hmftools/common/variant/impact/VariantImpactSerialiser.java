package com.hartwig.hmftools.common.variant.impact;

import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.fromSnpEffEnrichedVariant;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

// methods for reading from and writing to VCFs
public final class VariantImpactSerialiser
{
    public static final String VAR_IMPACT = "IMPACT";

    public static VCFHeader writeHeader(final VCFHeader header)
    {
        StringJoiner fields = new StringJoiner(", ");
        fields.add("Gene");
        fields.add("Transcript");
        fields.add("CanonicalEffect");
        fields.add("CanonicalCodingEffect");
        fields.add("SpliceRegion");
        fields.add("HgvsCodingImpact");
        fields.add("HgvsProteinImpact");
        fields.add("OtherReportableEffects");
        fields.add("WorstCodingEffect");
        fields.add("GenesAffected");

        header.addMetaDataLine(new VCFInfoHeaderLine(
                VAR_IMPACT, fields.length(), VCFHeaderLineType.String, String.format("PAVE Variant Impact [{}]", fields.toString())));
        return header;
    }

    public static void writeImpactDetails(
            final VariantContext context, final VariantImpact variantImpact)
    {
        context.getCommonInfo().putAttribute(VAR_IMPACT, VariantImpactSerialiser.toVcfData(variantImpact), true);
    }

    public static List<String> toVcfData(final VariantImpact impact)
    {
        return Lists.newArrayList(
                impact.CanonicalGeneName, impact.CanonicalTranscript, writeEffect(impact.CanonicalEffect),
                String.valueOf(impact.CanonicalCodingEffect), String.valueOf(impact.CanonicalSpliceRegion),
                impact.CanonicalHgvsCodingImpact, impact.CanonicalHgvsProteinImpact,
                impact.OtherReportableEffects, String.valueOf(impact.WorstCodingEffect), String.valueOf(impact.GenesAffected));
    }

    public static VariantImpact fromVariantContext(final VariantContext context)
    {
        if(context.hasAttribute(VAR_IMPACT))
            return fromAttributeValues(context.getAttributeAsStringList(VAR_IMPACT, ""));

        // revert to SnpEff until migration is complete
        return fromSnpEffEnrichedVariant(context);
    }

    public static VariantImpact fromAttributeValues(final List<String> impactValues)
    {
        int index = 0;
        String canonicalGeneName = impactValues.get(index++);
        String canonicalTranscript = impactValues.get(index++);
        String canonicalEffect = impactValues.get(index++);
        CodingEffect canonicalCodingEffect = CodingEffect.valueOf(impactValues.get(index++));

        boolean canonicalSpliceRegion = Boolean.parseBoolean(impactValues.get(index++));
        String canonicalHgvsCodingImpact = impactValues.get(index++);
        String canonicalHgvsProteinImpact = impactValues.get(index++);

        String otherReportableEffects = impactValues.get(index++);

        CodingEffect worstCodingEffect = CodingEffect.valueOf(impactValues.get(index++));
        int genesAffected = Integer.parseInt(impactValues.get(index++));

        return new VariantImpact(
                canonicalGeneName, canonicalTranscript, canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, canonicalSpliceRegion, otherReportableEffects, worstCodingEffect, genesAffected);

    }

    private static String writeEffect(final String effect)
    {
        return effect.replace("; ", "&").replace(" ", "_");
    }

}
