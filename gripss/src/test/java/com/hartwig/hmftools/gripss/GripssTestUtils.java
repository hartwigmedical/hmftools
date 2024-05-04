package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.region.ExcludedRegions.POLY_G_REGIONS_V37;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.BEID;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.BEIDL;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.CIRPOS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.EVENT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_AS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASRP;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASSR;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BAQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BSC;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BUM;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BUMQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_CAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_RAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_RPQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_SRQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.INDEL_COUNT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.PAR_ID;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.READ_PAIRS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SGL_FRAG_COUNT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.VcfIdGenerator.vcfId;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_SHORT_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_END;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_POINT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_TUMOR_AF_SGL;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_PON_DISTANCE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.PMS2_V37;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.gridss.GridssSvFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.filter.CompoundFilter;

public class GripssTestUtils
{
    public static final String LINE_INSERT_SEQ_A = "AAAAAAAAAAAAAAAAAAAA";
    public static final String LINE_INSERT_SEQ_T = "TTTTTTTTTTTTTTTTTTTT";

    public static final double DEFAULT_QUAL = 5000;

    public static GridssSvFactory defaultSvFactory()
    {
        return new GridssSvFactory(new CompoundFilter(false));
    }

    public static SvData createSv(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String insSeq, final GenotypeIds genotypeIds)
    {
        return createSv(eventId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, insSeq, genotypeIds,
                null, null, null);
    }

    public static SvData createSv(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String insSeq, final GenotypeIds genotypeIds, final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides,
            final Map<String,Object> tumorOverrides)
    {
        String ref = "A";

        VariantContext[] contexts = createSvBreakends(
                eventId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, ref, insSeq, commonOverrides, refOverrides, tumorOverrides);

        StructuralVariant sv = defaultSvFactory().createSV(contexts[SE_START], contexts[SE_END]);
        return new SvData(sv, genotypeIds);
    }

    public static SvData createSv(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String insSeq, final GenotypeIds genotypeIds, final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        String ref = "A";

        VariantContext[] contexts = createSvBreakends(
                eventId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, ref, insSeq, attributesStart, attributesEnd);

        StructuralVariant sv = defaultSvFactory().createSV(contexts[SE_START], contexts[SE_END]);
        return new SvData(sv, genotypeIds);
    }

    public static SvData createSgl(
            final String eventId, final String chromosome, int position, byte orientation, final String insSeq, final GenotypeIds genotypeIds)
    {
        return createSgl(
                eventId, chromosome, position, orientation, insSeq, genotypeIds, null, null, null);
    }

    public static SvData createSgl(
            final String eventId, final String chromosome, int position, byte orientation, final String insSeq, final GenotypeIds genotypeIds,
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        String ref = "A";

        VariantContext context = createSglBreakend(
                eventId, chromosome, position, orientation, ref, insSeq, commonOverrides, refOverrides, tumorOverrides);

        StructuralVariant sv = defaultSvFactory().createSingleBreakend(context);
        return new SvData(sv, genotypeIds);
    }

    public static VariantContext[] createSvBreakends(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String ref, final String insSeq)
    {
        return createSvBreakends(
                eventId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, ref, insSeq,
                null, null, null);
    }

    public static VariantContext[] createSvBreakends(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String ref, final String insSeq, final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides,
            final Map<String,Object> tumorOverrides)
    {
        String vcfStart = vcfId(eventId, true);
        String vcfEnd = vcfId(eventId, false);

        VariantContext[] pair = new VariantContext[SE_PAIR];

        String altStart = formPairedAltString(ref, insSeq, chrEnd, posEnd, orientStart, orientEnd);
        String altEnd = formPairedAltString(ref, insSeq, chrStart, posStart, orientEnd, orientStart);

        pair[SE_START] = createBreakend(vcfStart, chrStart, posStart, ref, altStart, vcfEnd, commonOverrides, refOverrides, tumorOverrides);
        pair[SE_END] = createBreakend(vcfEnd, chrEnd, posEnd, ref, altEnd, vcfStart, commonOverrides, refOverrides, tumorOverrides);

        return pair;
    }

    public static VariantContext[] createSvBreakends(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String ref, final String insSeq, final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        String vcfStart = vcfId(eventId, true);
        String vcfEnd = vcfId(eventId, false);

        VariantContext[] pair = new VariantContext[SE_PAIR];

        String altStart = formPairedAltString(ref, insSeq, chrEnd, posEnd, orientStart, orientEnd);
        String altEnd = formPairedAltString(ref, insSeq, chrStart, posStart, orientEnd, orientStart);

        pair[SE_START] = createBreakend(vcfStart, chrStart, posStart, ref, altStart, vcfEnd, attributesStart, null, null);
        pair[SE_END] = createBreakend(vcfEnd, chrEnd, posEnd, ref, altEnd, vcfStart, attributesEnd, null, null);

        return pair;
    }

    public static VariantContext createSglBreakend(
            final String eventId, final String chromosome, int position, byte orientation, final String ref, final String insSeq)
    {
        return createSglBreakend(eventId, chromosome, position, orientation, ref, insSeq, null, null, null);
    }

    public static VariantContext createSglBreakend(
            final String eventId, final String chromosome, int position, byte orientation, final String ref, final String insSeq,
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        String vcfId = vcfId(eventId, true);

        String alt = formSingleAltString(ref, insSeq, orientation);

        return createBreakend(vcfId, chromosome, position, ref, alt, null, commonOverrides, refOverrides, tumorOverrides);
    }

    public static VariantContext createBreakend(
            final String vcfId, final String chromosome, int position, final String ref, final String alt, final String mateId,
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        double qual = DEFAULT_QUAL;

        Map<String,Object> commonAttributes = makeCommonAttributes(vcfId, qual);

        if(commonOverrides != null)
            commonAttributes.putAll(commonOverrides);

        if(mateId != null)
            commonAttributes.put(PAR_ID, mateId);

        Map<String,Object> refAttributes = makeGenotypeAttributes(qual);
        Map<String,Object> tumorAttributes = makeGenotypeAttributes(qual);

        // defaults to indicate a somatic variant
        tumorAttributes.put(SV_FRAG_COUNT, 50);
        tumorAttributes.put(SGL_FRAG_COUNT, 50);
        tumorAttributes.put(GRIDSS_BSC, 100);
        tumorAttributes.put(SPLIT_READS, 1);

        if(refOverrides != null)
            refAttributes.putAll(refOverrides);

        if(tumorOverrides != null)
            tumorAttributes.putAll(tumorOverrides);

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

        double logError = -(qual / 10.0);

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

    public static Map<String,Object> makeCommonAttributes(final String vcfId, double qual)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        String eventId = vcfId.substring(0, vcfId.length() - 1);

        // quals
        attributes.put(QUAL, qual);
        attributes.put(GRIDSS_BQ, qual);
        attributes.put(GRIDSS_BAQ, qual);
        attributes.put(GRIDSS_SRQ, qual);
        attributes.put(GRIDSS_RPQ, qual);
        attributes.put(GRIDSS_BUMQ, 0);

        // read counts
        attributes.put(SPLIT_READS, 1);
        attributes.put(SV_FRAG_COUNT, 100);
        attributes.put(READ_PAIRS, 1);
        attributes.put(GRIDSS_ASRP, 1);
        attributes.put(GRIDSS_ASSR, 1);
        attributes.put(GRIDSS_BUM, 0);
        attributes.put(SGL_FRAG_COUNT, 100);

        attributes.put(INDEL_COUNT, 0);
        attributes.put(BEID, "");
        attributes.put(BEIDL, "");
        attributes.put(HOMSEQ, "");

        attributes.put(GRIDSS_AS, 0);
        attributes.put(GRIDSS_CAS, 0);
        attributes.put(GRIDSS_RAS, 0);

        attributes.put(EVENT, eventId);
        attributes.put(STRAND_BIAS, 0.5);
        attributes.put(REF_DEPTH_PAIR, 1);
        attributes.put(CIPOS, Lists.newArrayList(0, 0));
        attributes.put(CIRPOS, Lists.newArrayList(0, 0));
        attributes.put("SVTYPE", "BND");

        return attributes;
    }

    public static Map<String,Object> makeGenotypeAttributes(double qual)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        // quals
        attributes.put(QUAL, qual);
        attributes.put(GRIDSS_BQ, qual);
        attributes.put(GRIDSS_BAQ, qual);
        attributes.put(GRIDSS_SRQ, qual);
        attributes.put(GRIDSS_RPQ, qual);
        attributes.put(GRIDSS_BUMQ, 0);

        // read counts
        attributes.put(SPLIT_READS, 0);
        attributes.put(SV_FRAG_COUNT, 1);
        attributes.put(READ_PAIRS, 1);
        attributes.put(GRIDSS_ASRP, 1);
        attributes.put(GRIDSS_ASSR, 1);
        attributes.put(GRIDSS_BUM, 0);
        attributes.put(SGL_FRAG_COUNT, 1);
        attributes.put(GRIDSS_BSC, 1);
        attributes.put(REF_DEPTH, 10);

        // other
        attributes.put(INDEL_COUNT, 0);
        attributes.put(BEID, "");
        attributes.put(BEIDL, "");
        attributes.put(HOMSEQ, "");

        attributes.put(GRIDSS_AS, 0);
        attributes.put(GRIDSS_CAS, 0);
        attributes.put(GRIDSS_RAS, 0);

        attributes.put(REF_DEPTH_PAIR, 1);

        return attributes;
    }

    public static FilterConstants defaultFilterConstants()
    {
        return new FilterConstants(
                DEFAULT_HARD_MIN_TUMOR_QUAL,
                DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT,
                DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT,
                DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT,
                DEFAULT_MIN_NORMAL_COVERAGE,
                DEFAULT_MIN_TUMOR_AF_SGL,
                DEFAULT_MIN_TUMOR_AF,
                DEFAULT_MAX_SHORT_STRAND_BIAS,
                DEFAULT_MIN_QUAL_BREAK_END,
                DEFAULT_MIN_QUAL_BREAK_POINT,
                DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION,
                DEFAULT_MAX_HOM_LENGTH_SHORT_INV,
                DEFAULT_MIN_LENGTH,
                DEFAULT_PON_DISTANCE,
                POLY_G_REGIONS_V37,
                PMS2_V37,
                false,
                30,
                0.03,
                0.005);
    }

    public static Map<String,Object> buildLinkAttributes(final String beid, final String beidl)
    {
        Map<String,Object> attributes = Maps.newHashMap();
        attributes.put(GRIDSS_AS, 2); // set automatically from assembly strings

        String[] beids = beid.split(",");
        String[] beidls = beidl.split(",");
        attributes.put(BEID, beids);
        attributes.put(BEIDL, beidls);
        return attributes;
    }

    public static void loadSvDataCache(final SvDataCache dataCache, final List<SvData> svDataList)
    {
        dataCache.clear();
        svDataList.forEach(x -> dataCache.addSvData(x));
        dataCache.buildBreakendMap();
    }
}
