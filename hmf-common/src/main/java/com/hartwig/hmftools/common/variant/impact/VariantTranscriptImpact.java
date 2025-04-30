package com.hartwig.hmftools.common.variant.impact;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VariantTranscriptImpact
{
    public final String GeneId;
    public final String GeneName;
    public final String Transcript;
    public final String Effects;
    public final boolean SpliceRegion;
    public final String HgvsCoding;
    public final String HgvsProtein;
    public final String RefSeqId;
    public final int AffectedExon;
    public final int AffectedCodon;

    public VariantTranscriptImpact(
            final String geneId, final String geneName, final String transcript, final String effects,
            final boolean spliceRegion, final String hgvsCoding, final String hgvsProtein, final String refSeqId,
            final int affectedExon, final int affectedCodon)
    {
        GeneId = geneId;
        GeneName = geneName;
        Transcript = transcript;
        Effects = effects;
        SpliceRegion = spliceRegion;
        HgvsCoding = hgvsCoding;
        HgvsProtein = hgvsProtein;
        RefSeqId = refSeqId;
        AffectedExon = affectedExon;
        AffectedCodon = affectedCodon;
    }

    // serialisation
    public static final String VAR_TRANS_IMPACT_ANNOTATION = "PAVE_TI";

    // the in the VCF, transcript impacts are separated by ',', the components by '|' and the effects by '&"
    public static final String VAR_TRANS_IMPACT_DELIM = ",";
    public static final String VAR_TRANS_IMPACT_ITEM_DELIM = "|";

    public static void writeHeader(final VCFHeader header)
    {
        StringJoiner fields = new StringJoiner("|");
        List<String> fieldItems = Lists.newArrayList(
                "Gene", "GeneName", "Transcript", "Effects", "SpliceRegion", "HGVS.c", "HGVS.p", "RefSeqId", "Exon", "Codon");
        fieldItems.forEach(fields::add);

        header.addMetaDataLine(new VCFInfoHeaderLine(
                VAR_TRANS_IMPACT_ANNOTATION, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
                String.format("Transcript impact [%s]", fields)));
    }

    public static void writeVcfData(final VariantContext context, final List<VariantTranscriptImpact> transImpacts)
    {
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_DELIM);
        transImpacts.forEach(x -> sj.add(x.toVcfData()));
        context.getCommonInfo().putAttribute(VAR_TRANS_IMPACT_ANNOTATION, sj.toString(), true);
    }

    public static List<VariantTranscriptImpact> fromVariantContext(final VariantContext variant)
    {
        if(!variant.hasAttribute(VAR_TRANS_IMPACT_ANNOTATION))
            return Collections.EMPTY_LIST;

        String[] impactsStr = variant.getAttributeAsString(VAR_TRANS_IMPACT_ANNOTATION, "").split(VAR_TRANS_IMPACT_DELIM, -1);

        List<VariantTranscriptImpact> transImpacts = Lists.newArrayList();

        for(String impactStr : impactsStr)
        {
            transImpacts.add(fromVcfData(impactStr));
        }

        return transImpacts;
    }

    public static VariantTranscriptImpact fromVcfData(final String data)
    {
        String[] items = data.split("\\" + VAR_TRANS_IMPACT_ITEM_DELIM, -1);

        // RefSeqId and affect codon and exon added in v1.8
        String refSeqId = "";
        int affectedExon = 0;
        int affectedCodon = 0;

        if(items.length >= 10)
        {
            refSeqId = items[7];
            affectedExon = Integer.parseInt(items[8]);
            affectedCodon = Integer.parseInt(items[9]);
        }

        return new VariantTranscriptImpact(
                items[0], items[1], items[2], items[3], Boolean.parseBoolean(items[4]), items[5], items[6], refSeqId, affectedExon, affectedCodon);
    }

    private String toVcfData()
    {
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_ITEM_DELIM);
        sj.add(GeneId);
        sj.add(GeneName);
        sj.add(Transcript);
        sj.add(Effects);
        sj.add(String.valueOf(SpliceRegion));
        sj.add(HgvsCoding);
        sj.add(HgvsProtein);
        sj.add(RefSeqId);
        sj.add(String.valueOf(AffectedExon));
        sj.add(String.valueOf(AffectedCodon));
        return sj.toString();
    }

    public String toString() { return toVcfData(); }
}
