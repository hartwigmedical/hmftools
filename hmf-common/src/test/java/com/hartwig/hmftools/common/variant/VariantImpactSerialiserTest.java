package com.hartwig.hmftools.common.variant;

import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.junit.Test;

public class VariantImpactSerialiserTest
{
    @Test
    public void testSerialise()
    {
        VariantImpact variantImpact = new VariantImpact(
                "GENE", "TRANS_ID", "frameshift_variant&splice_region",
                CodingEffect.NONSENSE_OR_FRAMESHIFT, "HgvsCoding", "HgvsProtein",
                true, "", CodingEffect.MISSENSE, 1);

        List<String> impactValues = VariantImpactSerialiser.toVcfData(variantImpact);
        VariantImpact recreated = VariantImpactSerialiser.fromAttributeValues(impactValues);

        assertTrue(variantImpact.equals(recreated));
    }
}
