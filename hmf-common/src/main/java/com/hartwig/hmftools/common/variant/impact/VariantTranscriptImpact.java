package com.hartwig.hmftools.common.variant.impact;

import java.util.List;
import java.util.StringJoiner;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VariantTranscriptImpact
{
    public final String GeneId;
    public final String GeneName;
    public final String Transcript;
    public final String Effects;
    public final String HgvsCoding;
    public final String HgvsProtein;

    public VariantTranscriptImpact(
            final String geneId, final String geneName, final String transcript, final String effects,
            final String hgvsCoding, final String hgvsProtein)
    {
        GeneId = geneId;
        GeneName = geneName;
        Transcript = transcript;
        Effects = effects;
        HgvsCoding = hgvsCoding;
        HgvsProtein = hgvsProtein;
    }

    // serialisation
    public static final String VAR_TRANS_IMPACT_ANNOATATION = "VTI";
    public static final int VAR_TRANS_IMPACT_PARTS = 8;

    public static final String VAR_TRANS_IMPACT_DELIM = ",";
    public static final String VAR_TRANS_IMPACT_ITEM_DELIM = "|";
    public static final String VAR_TRANS_IMPACT_EFFECTS_DELIM = "&";

    public static void writeHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(
                VAR_TRANS_IMPACT_ANNOATATION, VAR_TRANS_IMPACT_PARTS, VCFHeaderLineType.String,
                "Transcript impact [Allele | Effects | GeneName | GeneId | Transcript | ExonRank | HGVS.c | HGVS.p]"));
    }

    public static void writeVcfData(
            final VariantContext context, final List<VariantTranscriptImpact> transImpacts)
    {
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_DELIM);
        transImpacts.forEach(x -> sj.add(x.toVcfData()));
        context.getCommonInfo().putAttribute(VAR_TRANS_IMPACT_ANNOATATION, sj.toString(), true);
    }

    public VariantTranscriptImpact fromVcfData(final String data)
    {
        String[] items = data.split(VAR_TRANS_IMPACT_ITEM_DELIM, -1);
        return new VariantTranscriptImpact(items[0], items[1], items[2], items[3], items[4], items[5]);
    }

    private String toVcfData()
    {
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_ITEM_DELIM);
        sj.add(GeneId);
        sj.add(GeneName);
        sj.add(Transcript);
        sj.add(Effects);
        sj.add(HgvsCoding);
        sj.add(HgvsProtein);
        return sj.toString();
    }

    public static String effectsToVcf(final List<String> effects)
    {
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_EFFECTS_DELIM);
        effects.forEach(x -> sj.add(x.toString()));
        return sj.toString();

    }
}
