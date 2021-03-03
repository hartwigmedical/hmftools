package com.hartwig.hmftools.common.drivercatalog.panel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineResourcesTest {

    @Test
    public void testBlacklistConsistency() throws IOException {
        testConsistency(GermlineResources.grch37Blacklist(), GermlineResources.grch38Blacklist());
    }

    @Test
    public void testWhiteConsistency() throws IOException {
        testConsistency(GermlineResources.grch37Whitelist(), GermlineResources.grch38Whitelist());
    }

    public void testConsistency(List<VariantContext> hg19list, List<VariantContext> hg38list )  {

        assertEquals(hg19list.size(), hg38list.size());

        for (int i = 0; i < hg19list.size(); i++) {
            VariantContext hg19 = hg19list.get(i);
            VariantContext hg38 = hg38list.get(i);

            assertEquals(hg19.getContig(), hg38.getContig().replace("chr", ""));
            assertEquals(hg19.getReference().getBaseString(), hg38.getReference().getBaseString());
            assertEquals(hg19.getAlternateAllele(0).getBaseString(), hg38.getAlternateAllele(0).getBaseString());
            assertTrue(hg38.getContig().contains("chr"));
        }

    }


}
