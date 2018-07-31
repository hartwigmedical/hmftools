package com.hartwig.hmftools.common.variant.cosmic;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.variant.VariantContextFromString;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class CosmicAnnotationFactoryTest {

    @Test
    public void canParseMultipleCosmicAnnotations() {
        final String line =
                "11\t62287501\tCOSM123;COSM456\tT\tC\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        VariantContext context = VariantContextFromString.decode(line);

        List<CosmicAnnotation> cosmicAnnotations = CosmicAnnotationFactory.fromContext(context);

        assertEquals(2, cosmicAnnotations.size());
        CosmicAnnotation annotation1 = cosmicAnnotations.get(0);
        assertEquals("GENE", annotation1.gene());
        assertEquals("TRANS1", annotation1.transcript());
        assertEquals("c.1A>G", annotation1.hgvsCoding());
        assertEquals("p.E1E", annotation1.hgvsProtein());
        assertEquals(1, annotation1.count());

        CosmicAnnotation annotation2 = cosmicAnnotations.get(1);
        assertEquals("GENE", annotation2.gene());
        assertEquals("TRANS2", annotation2.transcript());
        assertEquals("c.2A>G", annotation2.hgvsCoding());
        assertEquals("p.E2E", annotation2.hgvsProtein());
        assertEquals(1, annotation2.count());
    }
}