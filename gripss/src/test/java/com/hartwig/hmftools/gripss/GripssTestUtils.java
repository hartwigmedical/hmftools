package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_AS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BEID;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BEIDL;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BUM;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BUMQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BVF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CAS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIRPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_EVENT;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_IC;
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

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class GripssTestUtils
{
    public static final String CHR_1 = "1";
    public static final double DEFAULT_QUAL = 1000;

    public static VariantContext createBreakend(
            final String vcfId, final String chromosome, int position, final String ref, final String alt)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        double qual = DEFAULT_QUAL;

        Map<String,Object> attributes = Maps.newHashMap();

        String eventId = vcfId.substring(0, vcfId.length() - 1);

        attributes.put(VT_QUAL, qual);
        attributes.put(VT_SR, 1);
        attributes.put(VT_BQ, qual);
        attributes.put(VT_SRQ, qual);
        attributes.put(VT_VF, 1);
        attributes.put(VT_RP, 1);
        attributes.put(VT_IC, 0);
        attributes.put(VT_RPQ, qual);
        attributes.put(VT_REF, 1);
        attributes.put(VT_BEID, "");
        attributes.put(VT_BEIDL, "");
        attributes.put(VT_HOMSEQ, "");

        attributes.put(VT_AS, 0);
        attributes.put(VT_CAS, 0);
        attributes.put(VT_RAS, 0);

        attributes.put(VT_EVENT, eventId);
        attributes.put(VT_ASRP, 1);
        attributes.put(VT_SB, 0.5);
        attributes.put(VT_BVF, 1);
        attributes.put(VT_BUM, 0);
        attributes.put(VT_BUMQ, 0);
        attributes.put(VT_REFPAIR, 1);
        attributes.put(VT_CIPOS, 1);
        attributes.put(VT_CIRPOS, 1);
        attributes.put("SVTYPE", "BND");
        // attributes.put(VT_REALIGN, );

        Map<String,Object> refAttributes = Maps.newHashMap(attributes);
        Map<String,Object> tumorAttributes = Maps.newHashMap(attributes);

        GenotypeBuilder genotypeBuilder = new GenotypeBuilder();
        Genotype gtNormal = genotypeBuilder
                .attributes(refAttributes)
                .name(TEST_REF_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        Genotype gtTumor = genotypeBuilder
                .attributes(tumorAttributes)
                .name(TEST_SAMPLE_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        GenotypesContext genotypesContext = GenotypesContext.create(gtNormal, gtTumor);


        String filters = "";

        // String source, String ID, String contig, long start, long stop, Collection<Allele> alleles, GenotypesContext genotypes,
        // double log10PError, Set<String> filters, Map<String, Object> attributes,

        double logError = -(qual / 100.0)
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
}
