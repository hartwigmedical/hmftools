package com.hartwig.hmftools.sage.filters;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_RELATIVE_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.populateDefaultLevels;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_GERMLINE_MAX_AD;
import static com.hartwig.hmftools.sage.tinc.VariantCache.ChromosomeTask.checkFilters;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.tinc.FilterReason;
import com.hartwig.hmftools.sage.tinc.TincAnalyser;
import com.hartwig.hmftools.sage.tinc.TincCalculator;
import com.hartwig.hmftools.sage.tinc.VariantData;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class TincTest
{
    @Test
    public void testTincCalculations()
    {
        List<VariantData> variants = Lists.newArrayList();

        // too few variants
        double tincLevel = TincCalculator.calculate(variants, populateDefaultLevels());
        assertEquals(0, tincLevel, 0.001);

        // create different levels of ref support
        int refDepth = 30;
        int tumorDepth = 100;
        int tumorFrags = 30;

        for(int i = 0; i < 20; ++i)
        {
            for(int j = 0; j <= TINC_GERMLINE_MAX_AD; ++j)
            {
                variants.add(createTincVariant(refDepth, j, tumorDepth, tumorFrags));
            }
        }

        tincLevel = TincCalculator.calculate(variants, populateDefaultLevels());

        assertEquals(0.75, tincLevel, 0.001);
    }

    @Test
    public void testTincRecovery()
    {
        FilterConfig filterConfig = new FilterConfig();

        List<VariantData> variants = Lists.newArrayList();

        int refDepth = 30;
        int tumorDepth = 100;
        int tumorFrags = 30;

        VariantData var1 = createTincVariantFromContext(refDepth, 5, tumorDepth, tumorFrags, 30);
        var1.Context.getCommonInfo().addFilter(MAX_GERMLINE_RELATIVE_QUAL.filterName());
        var1.Context.getCommonInfo().addFilter(MAX_GERMLINE_VAF.filterName());
        var1.Context.getCommonInfo().addFilter(MAX_GERMLINE_ALT_SUPPORT.filterName());
        variants.add(var1);

        double tincLevel = 0.2;

        TincAnalyser.recoverVariants(filterConfig, variants, tincLevel);

        assertTrue(var1.newFilters().isEmpty());

        tincLevel = 0.001;

        TincAnalyser.recoverVariants(filterConfig, variants, tincLevel);

        assertEquals(2, var1.newFilters().size());
        assertTrue(var1.newFilters().contains(MAX_GERMLINE_RELATIVE_QUAL));
        assertTrue(var1.newFilters().contains(MAX_GERMLINE_VAF));
    }

    @Test
    public void testTincFilters()
    {
        double mapQualFactor = 30;
        VariantData variant = createTincVariantFromContext(30, 10, 100, 50, mapQualFactor);

        FilterReason filterReason = checkFilters(variant);

        assertEquals(FilterReason.GERMLINE_AF, filterReason);

        variant = createTincVariantFromContext(30, 1, 100, 50, mapQualFactor);
        variant.setPonFiltered();
        filterReason = checkFilters(variant);
        assertEquals(FilterReason.PON, filterReason);

        variant = createTincVariantFromContext(30, 1, 100, 50, mapQualFactor);
        variant.setGnomadFrequency(0.01);
        filterReason = checkFilters(variant);
        assertEquals(FilterReason.GNOMAD, filterReason);

        variant = createTincVariantFromContext(30, 1, 100, 50, mapQualFactor);
        variant.Context.getGenotype(0).getExtendedAttributes().put(AVG_RECALIBRATED_BASE_QUAL, "25,25");
        filterReason = checkFilters(variant);
        assertEquals(FilterReason.LOW_ABQ, filterReason);

        variant = createTincVariantFromContext(30, 1, 100, 50, 10);
        filterReason = checkFilters(variant);
        assertEquals(FilterReason.MQF, filterReason);
    }

    private static VariantData createTincVariant(
            final int referenceDepth, final int referenceAltFrags, final int tumorDepth, final int tumorAltFrags)
    {
        return new VariantData(CHR_1, 100, "A", "C", referenceDepth, referenceAltFrags, tumorDepth, tumorAltFrags);
    }

    private static final String TUMOR_SAMPLE = "TUMOR_ID";
    private static final String REF_SAMPLE = "REF_ID";

    private static VariantData createTincVariantFromContext(
            final int referenceDepth, final int referenceAltFrags, final int tumorDepth, final int tumorAltFrags, double mapQualFactor)
    {
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(REF_SAMPLE);
        genotypeBuilder.AD(new int[] { 0, referenceAltFrags });
        genotypeBuilder.DP(referenceDepth);
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        genotypeBuilder.attribute(AVG_RECALIBRATED_BASE_QUAL, "30,30");

        int refFrags = referenceDepth - referenceAltFrags;
        genotypeBuilder.attribute(READ_CONTEXT_COUNT, format("%d,0,0,0,%d,%d", referenceAltFrags, refFrags, referenceDepth));
        genotypeBuilder.attribute(READ_CONTEXT_QUALITY, format("%d,0,0,0,%d,%d", referenceAltFrags, refFrags, referenceDepth));

        Genotype refGenotype = genotypeBuilder.make();

        genotypeBuilder = new GenotypeBuilder(TUMOR_SAMPLE);
        genotypeBuilder.AD(new int[] { referenceAltFrags, tumorAltFrags });
        genotypeBuilder.DP(tumorDepth);
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        genotypeBuilder.attribute(AVG_RECALIBRATED_BASE_QUAL, "30,30");

        refFrags = tumorDepth - tumorAltFrags;
        genotypeBuilder.attribute(READ_CONTEXT_COUNT, format("%d,0,0,0,%d,%d", tumorAltFrags, refFrags, tumorDepth));
        genotypeBuilder.attribute(READ_CONTEXT_QUALITY, format("%d,0,0,0,%d,%d", tumorAltFrags, refFrags, tumorDepth));

        Genotype tumorGenotype = genotypeBuilder.make();

        List<Allele> alleles = List.of(Allele.create("A", true), Allele.create("C", false));

        VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(CHR_1).start(100).computeEndFromAlleles(alleles, 100);
        builder.alleles(alleles);

        builder.genotypes(List.of(refGenotype, tumorGenotype));

        builder.attribute(MAP_QUAL_FACTOR, mapQualFactor);

        VariantContext variantContext = builder.make(true);

        GenotypeIds genotypeIds = new GenotypeIds(0, 1, REF_SAMPLE, TUMOR_SAMPLE);

        return new VariantData(variantContext, genotypeIds);
    }
}
