package com.hartwig.hmftools.sigs;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SmallVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.sigs.fitter.SigVariant;
import com.hartwig.hmftools.sigs.fitter.VariantFilters;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.HotspotType.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class SignaturesTest
{
    @Test
    public void testVariantFilters()
    {
        SimpleVariant simpleVariant = new SimpleVariant(CHR_1, 100, "A", "T");
        VariantContext variantContext = buildVariantContext(simpleVariant);

        SmallVariantFactory variantFactory = SmallVariantFactory.passOnlyInstance();

        SmallVariant smallVariant = variantFactory.createVariant(TEST_SAMPLE_ID, variantContext).get();
        SigVariant variant = new SigVariant(smallVariant, variantContext);

        VariantFilters filters = new VariantFilters(null, null, null, null, true);

        variantContext.getCommonInfo().putAttribute(PURPLE_AF, 0.04);

        assertFalse(filters.passesFilters(variant));

        variantContext.getCommonInfo().removeAttribute(PURPLE_AF);
        variantContext.getCommonInfo().putAttribute(PURPLE_AF, 0.1);

        assertTrue(filters.passesFilters(variant));

        variantContext.getCommonInfo().putAttribute(GNOMAD_FREQ, 0.001);

        assertFalse(filters.passesFilters(variant));

        variantContext.getCommonInfo().removeAttribute(GNOMAD_FREQ);

        assertTrue(filters.passesFilters(variant));

        variantContext.getCommonInfo().putAttribute(HOTSPOT_FLAG, true);

        smallVariant = variantFactory.createVariant(TEST_SAMPLE_ID, variantContext).get();
        variant = new SigVariant(smallVariant, variantContext);

        assertFalse(filters.passesFilters(variant));

        variantContext.getCommonInfo().removeAttribute(HOTSPOT_FLAG);
        variantContext.getCommonInfo().putAttribute(PANEL_SOMATIC_LIKELIHOOD, SomaticLikelihood.LOW.toString());

        smallVariant = variantFactory.createVariant(TEST_SAMPLE_ID, variantContext).get();
        variant = new SigVariant(smallVariant, variantContext);

        assertFalse(filters.passesFilters(variant));
    }

    private static final String TEST_SAMPLE_ID = "Sample_001";

    private static VariantContext buildVariantContext(final SimpleVariant variant)
    {
        GenotypeBuilder genotypeBuilder = new GenotypeBuilder(TEST_SAMPLE_ID);
        genotypeBuilder.AD(new int[] { 10, 10 });
        genotypeBuilder.DP(100);
        genotypeBuilder.alleles(Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL));
        Genotype genotype = genotypeBuilder.make();

        List<Allele> alleles = List.of(Allele.create(variant.Ref, true), Allele.create(variant.Alt, false));

        VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(variant.Chromosome).start(variant.Position);
        builder.alleles(alleles);
        builder.computeEndFromAlleles(alleles, variant.Position);

        builder.genotypes(List.of(genotype));

        return builder.make();
    }

}
