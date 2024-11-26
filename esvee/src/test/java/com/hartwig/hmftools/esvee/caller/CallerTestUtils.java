package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.filter.CompoundFilter;

public class CallerTestUtils
{
    public static final String TEST_REF_ID = "REF_ID";
    public static final String TEST_SAMPLE_ID = "SAMPLE_ID";
    public static final int DEFAULT_QUAL = 60;

    private static final StructuralVariantFactory SV_FACTOR = new StructuralVariantFactory(new CompoundFilter(false));

    public static GenotypeIds createGenotypeIds()
    {
        return new GenotypeIds(0, 10, TEST_REF_ID, TEST_SAMPLE_ID);
    }

    public static Variant createSv(
            final String vcfId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String insSeq)
    {
        return createSv(
                vcfId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, insSeq,
                Collections.emptyMap(), Collections.emptyMap(), Collections.emptyMap());
    }

    public static Variant createSv(
            final String vcfId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String insSeq, final Map<String,Object> attributes, final Map<String,Object> refAttributes, final Map<String,Object> tumorAttributes)
    {
        String ref = "A";

        String altStart = formPairedAltString(ref, insSeq, chrEnd, posEnd, orientStart, orientEnd);
        String altEnd = formPairedAltString(ref, insSeq, chrStart, posStart, orientEnd, orientStart);

        String vcfIdStart = vcfId + "_start";
        String vcfIdEnd = vcfId + "_end";

        VariantContext contextStart = createContext(
                vcfIdStart, chrStart, posStart, ref, altStart, vcfIdEnd, attributes, refAttributes, tumorAttributes);

        GenotypeIds genotypeIds = createGenotypeIds();

        if(chrEnd == null)
        {
            StructuralVariant sv = SV_FACTOR.createSingleBreakend(contextStart);
            return new Variant(sv, genotypeIds);
        }

        VariantContext contextEnd = createContext(
                vcfIdEnd, chrEnd, posEnd, ref, altEnd, vcfIdStart, attributes, refAttributes, tumorAttributes);

        StructuralVariant sv = SV_FACTOR.createSV(contextStart, contextEnd);
        return new Variant(sv, genotypeIds);
    }

    public static VariantContext createContext(
            final String vcfId, final String chromosome, int position, final String ref, final String alt, final String mateId,
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        double defaultQual = DEFAULT_QUAL;

        Map<String,Object> commonAttributes = commonOverrides != null ? commonOverrides : makeCommonAttributes(vcfId, 50, 50, 50);

        //if(mateId != null && !commonAttributes.isEmpty())
        //    commonAttributes.put(PAR_ID, mateId);

        Map<String,Object> refAttributes = refOverrides != null ? refOverrides : makeGenotypeAttributes(0, 0, 30);
        Map<String,Object> tumorAttributes = tumorOverrides != null ? tumorOverrides : makeGenotypeAttributes(30, 20, 100);

        Genotype gtNormal = new GenotypeBuilder()
                .attributes(refAttributes)
                .name(TEST_REF_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        Genotype gtTumor = new GenotypeBuilder()
                .attributes(tumorAttributes)
                .name(TEST_SAMPLE_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        GenotypesContext genotypesContext = GenotypesContext.create(gtNormal, gtTumor);

        double logError = -(defaultQual / 10.0);

        return builder
                .source("SOURCE")
                .id(vcfId)
                .chr(chromosome)
                .start(position)
                .stop(position)
                .alleles(alleles)
                .genotypes(genotypesContext)
                .attributes(commonAttributes)
                .log10PError(logError)
                .unfiltered()
                .make(true);
    }

    public static Map<String,Object> makeCommonAttributes(final String vcfId, double qual, int splitFrags, int discFrags)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        attributes.put(QUAL, qual);
        attributes.put(SPLIT_FRAGS, splitFrags);
        attributes.put(DISC_FRAGS, discFrags);
        attributes.put(TOTAL_FRAGS, splitFrags + discFrags);
        attributes.put(AVG_FRAG_LENGTH, 500);

        return attributes;
    }

    public static Map<String,Object> makeGenotypeAttributes(int splitFrags, int discFrags, int refFrags)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        attributes.put(SPLIT_FRAGS, splitFrags);
        attributes.put(DISC_FRAGS, discFrags);
        attributes.put(TOTAL_FRAGS, splitFrags + discFrags);
        attributes.put(REF_DEPTH, refFrags);

        return attributes;
    }
}
