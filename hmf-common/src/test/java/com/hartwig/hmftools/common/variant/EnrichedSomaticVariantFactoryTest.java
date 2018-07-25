package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class EnrichedSomaticVariantFactoryTest {

    @Test
    public void revertToNormalCosmicIDIfCanonicalIsNoMatch() {
        final String line =
                "14\t105246551\trs121434592;COSM33765\tC\tT\t.\tPASS\tCOSM2ENST=COSM33765|AKT1_ENST00000349310|c.49G>A|p.E17K|519\tGT:AD:DP\t0/1:120,71:204";

        final String sample = "sample";
        final VariantContext context = VariantContextFromString.decode(sample, line);
        final Optional<SomaticVariant> optVariant = SomaticVariantFactory.passOnlyInstance().createVariant(sample, context);
        assert optVariant.isPresent();
        final SomaticVariant variant = optVariant.get();

        final ImmutableEnrichedSomaticVariant.Builder builder = SomaticVariantTestBuilderFactory.createEnriched();
        final Map<String, String> geneToTranscriptMap = Maps.newHashMap();
        geneToTranscriptMap.put("AKT1", "ENST0123456");
        final TranscriptAnnotationSelector selector = new TranscriptAnnotationSelector(geneToTranscriptMap);

        EnrichedSomaticVariantFactory.addCanonicalCosmicID(builder, variant, selector);

        final EnrichedSomaticVariant enrichedSomaticVariant = builder.build();

        assertEquals("COSM33765", enrichedSomaticVariant.canonicalCosmicID());
    }
}
