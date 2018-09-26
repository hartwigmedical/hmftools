package com.hartwig.hmftools.common.variant.filter;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.VariantContextFromString;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class NearIndelPonFilterTest {

    @Test
    public void testFindPonForwards() {
        final List<VariantContext> victims = Lists.newArrayList(create(87, "TAG", "A"),
                create(87, "TAGG", "A"),
                create(88, "TAG", "A"),
                create(100, "TAG", "A"),
                createPonFiltered(100, "TAG", "A"));

        assertFalse(NearIndelPonFilter.isIndelNearPon(0, victims));
        assertTrue(NearIndelPonFilter.isIndelNearPon(1, victims));
        assertTrue(NearIndelPonFilter.isIndelNearPon(2, victims));
        assertTrue(NearIndelPonFilter.isIndelNearPon(3, victims));
        assertFalse(NearIndelPonFilter.isIndelNearPon(4, victims));
    }

    @Test
    public void testFindPonBackwards() {
        final List<VariantContext> victims = Lists.newArrayList(createPonFiltered(100, "TAG", "A"),
                create(100, "TAG", "A"),
                create(112, "TAG", "A"),
                create(113, "TAG", "A"),
                create(113, "TAGGG", "A"));

        assertFalse(NearIndelPonFilter.isIndelNearPon(0, victims));
        assertTrue(NearIndelPonFilter.isIndelNearPon(1, victims));
        assertTrue(NearIndelPonFilter.isIndelNearPon(2, victims));
        assertFalse(NearIndelPonFilter.isIndelNearPon(3, victims));
        assertFalse(NearIndelPonFilter.isIndelNearPon(4, victims));
    }

    @Test
    public void testIgnoreSNP() {
        final List<VariantContext> victims =
                Lists.newArrayList(create(99, "TT", "A"), create(100, "T", "A"), createPonFiltered(100, "TAG", "A"));

        assertTrue(NearIndelPonFilter.isIndelNearPon(0, victims));
        assertFalse(NearIndelPonFilter.isIndelNearPon(1, victims));
    }

    @Test
    public void testNoPon() {
        final List<VariantContext> victims = Lists.newArrayList(create(99, "TT", "A"), create(100, "TAG", "A"));

        assertFalse(NearIndelPonFilter.isIndelNearPon(0, victims));
        assertFalse(NearIndelPonFilter.isIndelNearPon(1, victims));
    }

    @NotNull
    private static VariantContext create(int start, @NotNull final String ref, @NotNull final String alt) {
        final String line = "1\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tC\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }

    @NotNull
    private static VariantContext createPonFiltered(int start, @NotNull final String ref, @NotNull final String alt) {
        final String line = "1\t" + start + "\t.\t" + ref + "\t" + alt + "\t.\tGERMLINE_PON\t;SOMATIC_PON_HET=34;SOMATIC_PON_HOM=0\t";
        return VariantContextFromString.decode(line);
    }
}
