package com.hartwig.hmftools.common.variant.impact;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
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

    // in the VCF, the components of the variant impact are separated by ',' and effects are separated by '&'
    // other reportable effects are separated by '-' and their sub-details by '|'

    public static final String VAR_IMPACT_OTHER_REPORT_ITEM_DELIM = "|";
    public static final String VAR_IMPACT_OTHER_REPORT_DELIM = "-";

    public static VCFHeader writeHeader(final VCFHeader header)
    {
        StringJoiner fields = new StringJoiner(", ");

        List<String> fieldItems = Lists.newArrayList(
                "Gene", "Transcript", "CanonicalEffect", "CanonicalCodingEffect", "SpliceRegion", "HgvsCodingImpact",
        "HgvsProteinImpact", "OtherReportableEffects", "WorstCodingEffect", "GenesAffected");
        fieldItems.forEach(x -> fields.add(x));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                VAR_IMPACT, fieldItems.size(), VCFHeaderLineType.String, String.format("Variant Impact [%s]", fields.toString())));
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
                impact.CanonicalGeneName,
                impact.CanonicalTranscript,
                impact.CanonicalEffect,
                String.valueOf(impact.CanonicalCodingEffect),
                String.valueOf(impact.CanonicalSpliceRegion),
                impact.CanonicalHgvsCoding,
                impact.CanonicalHgvsProtein,
                impact.OtherReportableEffects,
                String.valueOf(impact.WorstCodingEffect),
                String.valueOf(impact.GenesAffected));
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
        if(impactValues.size() != 10)
        {
            return new VariantImpact(
                    "", "", "", UNDEFINED, "", "",
                    false, "", UNDEFINED, 0);
        }

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

    public static String toOtherReportableTransInfo(
            final String transName, final String hgvsCoding, final String hgvsProtein, final String effects, final CodingEffect codingEffect)
    {
        // eg ENST00000579755|c.209_210delCCinsTT|p.Pro70Leu|missense_variant|MISSENSE;
        StringJoiner sj = new StringJoiner(VAR_IMPACT_OTHER_REPORT_ITEM_DELIM);
        sj.add(transName);
        sj.add(hgvsCoding);
        sj.add(hgvsProtein);
        sj.add(effects);
        sj.add(codingEffect.toString());
        return sj.toString();
    }
}
