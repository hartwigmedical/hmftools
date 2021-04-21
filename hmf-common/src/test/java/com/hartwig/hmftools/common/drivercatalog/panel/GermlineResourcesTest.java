package com.hartwig.hmftools.common.drivercatalog.panel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineResourcesTest {

    @Test
    public void testBlacklistConsistency() throws IOException {
        testConsistency(GermlineResources.blacklist37(), GermlineResources.blacklist38());
    }

    @Test
    public void testWhitelistConsistency() throws IOException {
        testConsistency(GermlineResources.whitelist37(), GermlineResources.whitelist38());
    }

    private static void testConsistency(@NotNull List<VariantContext> v37list, @NotNull List<VariantContext> v38list )  {
        assertEquals(v37list.size(), v38list.size());

        for (int i = 0; i < v37list.size(); i++) {
            VariantContext v37 = v37list.get(i);
            VariantContext v38 = v38list.get(i);

            if (v38.getStart() == 99672916) {
                // rs776746 is defined such that it indicates a C>T snp in v37 and a T>C snp in v38
                assertEquals(v37.getReference().getBaseString(), v38.getAlternateAllele(0).getBaseString());
                assertEquals(v37.getAlternateAllele(0).getBaseString(), v38.getReference().getBaseString());
            } else {
                assertEquals(v37.getReference().getBaseString(), v38.getReference().getBaseString());
                assertEquals(v37.getAlternateAllele(0).getBaseString(), v38.getAlternateAllele(0).getBaseString());
            }

            assertEquals(v37.getContig(), v38.getContig().replace("chr", ""));
            assertTrue(v38.getContig().contains("chr"));
        }
    }
}
