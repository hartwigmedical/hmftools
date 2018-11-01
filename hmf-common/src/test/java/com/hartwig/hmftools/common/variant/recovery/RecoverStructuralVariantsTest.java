package com.hartwig.hmftools.common.variant.recovery;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class RecoverStructuralVariantsTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testMate() {
        assertEquals("17:59493156", RecoverStructuralVariants.mateLocation("C[17:59493156["));
        assertEquals("17:59493156", RecoverStructuralVariants.mateLocation("]17:59493156]C"));
    }

    @Test
    public void testClosestRegion() {
        final GenomeRegion start = GenomeRegionFactory.create("1", 1, 10000);
        final GenomeRegion middle1 = GenomeRegionFactory.create("1", 10001, 20000);
        final GenomeRegion middle2 = GenomeRegionFactory.create("1", 20001, 30000);
        final GenomeRegion end = GenomeRegionFactory.create("1", 30001, 40000);
        final List<GenomeRegion> regions = Lists.newArrayList(start, middle1, middle2, end);

        assertEquals(start, RecoverStructuralVariants.closest(1, regions));
        assertEquals(start, RecoverStructuralVariants.closest(5001, regions));
        assertEquals(middle1, RecoverStructuralVariants.closest(5002, regions));
        assertEquals(middle1, RecoverStructuralVariants.closest(15001, regions));
        assertEquals(middle2, RecoverStructuralVariants.closest(15002, regions));
        assertEquals(middle2, RecoverStructuralVariants.closest(25001, regions));
        assertEquals(end, RecoverStructuralVariants.closest(25002, regions));
        assertEquals(end, RecoverStructuralVariants.closest(40000, regions));
    }

    @Test
    public void testFilterFilter() {
        final String eligible =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAA	343.95	BPI.Filter.SRSupportZero;qual	SVTYPE=BND	GT\t./.";
        final String passing =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAA	343.95	.	SVTYPE=BND	GT\t./.";
        final String ponFiltered =
                "16	33040203	gridss81_16816o	G	GGTAAGAATCCGC[16:33040204[	343.95	PON	SVTYPE=BND	GT\t./.";
        final String afFiltered =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAA	343.95	af	SVTYPE=BND	GT\t./.";

        assertTrue(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(eligible)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(passing)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(ponFiltered)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(afFiltered)));
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

}
