package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIG;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.pave.annotation.ClinvarChrCache;
import com.hartwig.hmftools.common.perf.StringCache;

import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class ResourcesTest
{
    @Test
    public void testClinvar()
    {
        VariantData var1 = createVariant(CHR_1, 100, "A", "G");
        VariantData var2 = createVariant(CHR_1, 100, "A", "C");
        VariantData var3 = createVariant(CHR_1, 100, "A", "T");

        StringCache strCache = new StringCache();
        ClinvarChrCache clinvarCache = new ClinvarChrCache(CHR_1, strCache);

        String matchedSig = "pathogenic";
        String conflictStr = "conflict";
        clinvarCache.addEntry(99, "A", "T", matchedSig, conflictStr);
        clinvarCache.addEntry(100, "A", "AT", matchedSig, conflictStr);
        clinvarCache.addEntry(100, "A", "C", matchedSig, conflictStr);
        clinvarCache.addEntry(100, "A", "G", matchedSig, conflictStr);
        clinvarCache.addEntry(101, "A", "T", matchedSig, conflictStr);

        clinvarCache.annotateVariant(var1);
        clinvarCache.annotateVariant(var2);
        clinvarCache.annotateVariant(var3);

        assertTrue(var1.context().hasAttribute(CLNSIG));
        assertTrue(var2.context().hasAttribute(CLNSIG));
        assertFalse(var3.context().hasAttribute(CLNSIG));
    }

    public static VariantData createVariant( final String chromosome, int position, final String ref, final String alt)
    {
        VariantContext context = buildContext(chromosome, position, ref, alt);
        VariantData variant = new VariantData(chromosome, position, ref, alt);
        variant.setContext(context);
        return variant;
    }

    public static VariantContext buildContext( final String chromosome, int position, final String ref, final String alt)
    {
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create(ref, true));
        alleles.add(Allele.create(alt, false));

        Map<String,Object> commonAttributes = Maps.newHashMap();

        return builder
                .source("SOURCE")
                .chr(chromosome)
                .start(position)
                .stop(position)
                .alleles(alleles)
                // .genotypes(genotypesContext)
                .attributes(commonAttributes)
                // .log10PError(logError)
                .unfiltered()
                .make(true);
    }
}
