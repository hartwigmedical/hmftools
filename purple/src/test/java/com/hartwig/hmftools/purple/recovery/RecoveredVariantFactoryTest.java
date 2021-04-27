package com.hartwig.hmftools.purple.recovery;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class RecoveredVariantFactoryTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testMateLocation() {
        assertEquals("17:59493156", RecoveredVariantFactory.mateLocation("C[17:59493156["));
        assertEquals("17:59493156", RecoveredVariantFactory.mateLocation("]17:59493156]C"));
    }

    @Test
    public void testFilterFilter() {
        final String eligible =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAA	343.95	BPI.Filter.SRSupportZero;someFilter	SVTYPE=BND	GT\t./.";
        final String passing =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAA	343.95	.	SVTYPE=BND	GT\t./.";
        final String ponFiltered =
                "16	33040203	gridss81_16816o	G	GGTAAGAATCCGC[16:33040204[	343.95	PON	SVTYPE=BND	GT\t./.";

        assertTrue(RecoveredVariantFactory.isAppropriatelyFiltered(codec.decode(eligible)));
        assertFalse(RecoveredVariantFactory.isAppropriatelyFiltered(codec.decode(passing)));

        assertTrue(RecoveredVariantFactory.isAppropriatelyFiltered(codec.decode(ponFiltered)));
        for (String filter : RecoveredVariantFactory.DO_NOT_RESCUE) {
            final String ineligible = ponFiltered.replace("PON", filter);
            assertFalse(RecoveredVariantFactory.isAppropriatelyFiltered(codec.decode(ineligible)));
        }
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
