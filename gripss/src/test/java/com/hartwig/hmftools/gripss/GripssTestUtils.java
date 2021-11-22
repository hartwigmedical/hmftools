package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.VcfIdGenerator.vcfId;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_AS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_ASSR;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BEID;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BEIDL;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BSC;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BUM;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BUMQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BVF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CAS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIRPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_EVENT;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_IC;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_PAR_ID;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_RAS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_RP;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_RPQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_SB;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_SR;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_SRQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_VF;
import static com.hartwig.hmftools.gripss.common.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.gripss.common.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_SHORT_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_END;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_POINT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_PON_DISTANCE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.LINC_00486_V37;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.filters.FilterConstants;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class GripssTestUtils
{
    public static final String CHR_1 = "1";
    public static final String CHR_2 = "2";

    public static final double DEFAULT_QUAL = 1000;

    public static VariantContext[] createSvBreakends(
            final String eventId, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String ref, final String insSeq)
    {
        String vcfStart = vcfId(eventId, true);
        String vcfEnd = vcfId(eventId, false);

        VariantContext[] pair = new VariantContext[SE_PAIR];

        String altStart = formPairedAltString(ref, insSeq, chrEnd, posEnd, orientStart, orientEnd);
        String altEnd = formPairedAltString(ref, insSeq, chrStart, posStart, orientEnd, orientStart);

        pair[SE_START] = createBreakend(vcfStart, chrStart, posStart, ref, altStart, vcfEnd);
        pair[SE_END] = createBreakend(vcfEnd, chrEnd, posEnd, ref, altEnd, vcfStart);

        return pair;
    }

    public static VariantContext createSglBreakend(
            final String eventId, final String chromosome, int position, byte orientation, final String ref, final String insSeq)
    {
        String vcfId = vcfId(eventId, true);

        String alt = formSingleAltString(ref, insSeq, orientation);

        return createBreakend(vcfId, chromosome, position, ref, alt, null);
    }

    public static VariantContext createBreakend(
            final String vcfId, final String chromosome, int position, final String ref, final String alt, final String mateId)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        double qual = DEFAULT_QUAL;

        Map<String,Object> attributes = makeCommonAttributes(vcfId, qual);

        if(mateId != null)
            attributes.put(VT_PAR_ID, mateId);

        Map<String,Object> refAttributes = makeGenotypeAttributes(qual);
        Map<String,Object> tumorAttributes = makeGenotypeAttributes(qual);

        tumorAttributes.put(VT_VF, 100); // indicating a somatic variant
        tumorAttributes.put(VT_BVF, 100);
        tumorAttributes.put(VT_BSC, 100);

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

        String filters = "";

        double logError = -(qual / 10.0)
;
        return builder
                .source("SOURCE")
                .id(vcfId)
                .chr(chromosome)
                .start(position)
                .stop(position)
                .alleles(alleles)
                .genotypes(genotypesContext)
                .attributes(attributes)
                .log10PError(logError)
                .filter(filters)
                .make(true);
    }

    public static Map<String,Object> makeCommonAttributes(final String vcfId, double qual)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        String eventId = vcfId.substring(0, vcfId.length() - 1);

        // quals
        attributes.put(VT_QUAL, qual);
        attributes.put(VT_BQ, qual);
        attributes.put(VT_SRQ, qual);
        attributes.put(VT_RPQ, qual);
        attributes.put(VT_BUMQ, 0);

        // read counts
        attributes.put(VT_SR, 1);
        attributes.put(VT_VF, 100);
        attributes.put(VT_RP, 1);
        attributes.put(VT_ASRP, 1);
        attributes.put(VT_ASSR, 1);
        attributes.put(VT_BUM, 0);
        attributes.put(VT_BVF, 100);
        attributes.put(VT_REF, 1);

        attributes.put(VT_IC, 0);
        attributes.put(VT_BEID, "");
        attributes.put(VT_BEIDL, "");
        attributes.put(VT_HOMSEQ, "");

        attributes.put(VT_AS, 0);
        attributes.put(VT_CAS, 0);
        attributes.put(VT_RAS, 0);

        attributes.put(VT_EVENT, eventId);
        attributes.put(VT_SB, 0.5);
        attributes.put(VT_REFPAIR, 1);
        attributes.put(VT_CIPOS, Lists.newArrayList(0, 0));
        attributes.put(VT_CIRPOS, Lists.newArrayList(0, 0));
        attributes.put("SVTYPE", "BND");

        // attributes.put(VT_REALIGN, );

        return attributes;
    }

    public static Map<String,Object> makeGenotypeAttributes(double qual)
    {
        Map<String,Object> attributes = Maps.newHashMap();

        // quals
        attributes.put(VT_QUAL, qual);
        attributes.put(VT_BQ, qual);
        attributes.put(VT_SRQ, qual);
        attributes.put(VT_RPQ, qual);
        attributes.put(VT_BUMQ, 0);

        // read counts
        attributes.put(VT_SR, 1);
        attributes.put(VT_VF, 1);
        attributes.put(VT_RP, 1);
        attributes.put(VT_ASRP, 1);
        attributes.put(VT_ASSR, 1);
        attributes.put(VT_BUM, 0);
        attributes.put(VT_BVF, 1);
        attributes.put(VT_BSC, 1);

        // other
        attributes.put(VT_IC, 0);
        attributes.put(VT_REF, 1);
        attributes.put(VT_BEID, "");
        attributes.put(VT_BEIDL, "");
        attributes.put(VT_HOMSEQ, "");

        attributes.put(VT_AS, 0);
        attributes.put(VT_CAS, 0);
        attributes.put(VT_RAS, 0);

        attributes.put(VT_REFPAIR, 1);

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
                DEFAULT_MIN_TUMOR_AF,
                DEFAULT_MAX_SHORT_STRAND_BIAS,
                DEFAULT_MIN_QUAL_BREAK_END,
                DEFAULT_MIN_QUAL_BREAK_POINT,
                DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION,
                DEFAULT_MAX_HOM_LENGTH_SHORT_INV,
                DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL,
                DEFAULT_MIN_LENGTH,
                DEFAULT_PON_DISTANCE,
                LINC_00486_V37);
    }
}
