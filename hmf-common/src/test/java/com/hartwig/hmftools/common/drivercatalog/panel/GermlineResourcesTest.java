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
        testConsistency(GermlineResources.blacklist37(), GermlineResources.blacklist38());
    }

    @Test
    public void testWhiteConsistency() throws IOException {
        testConsistency(GermlineResources.whitelist37(), GermlineResources.whitelist38());
    }

    public void testConsistency(List<VariantContext> hg19list, List<VariantContext> hg38list )  {

        assertEquals(hg19list.size(), hg38list.size());

        for (int i = 0; i < hg19list.size(); i++) {
            VariantContext hg19 = hg19list.get(i);
            VariantContext hg38 = hg38list.get(i);

            if (hg38.getStart() == 99672916) {
                // rs776746 is defined such that it indicates a C>T snp in hg19 and a T>C snp in hg38
                assertEquals(hg19.getReference().getBaseString(), hg38.getAlternateAllele(0).getBaseString());
                assertEquals(hg19.getAlternateAllele(0).getBaseString(), hg38.getReference().getBaseString());
            } else {
                assertEquals(hg19.getReference().getBaseString(), hg38.getReference().getBaseString());
                assertEquals(hg19.getAlternateAllele(0).getBaseString(), hg38.getAlternateAllele(0).getBaseString());
            }


            assertEquals(hg19.getContig(), hg38.getContig().replace("chr", ""));
            assertTrue(hg38.getContig().contains("chr"));
        }

    }


}
