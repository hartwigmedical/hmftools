package com.hartwig.hmftools.common.variant.snpeff;

import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_CONSEQUENCE;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser.SNPEFF_IDENTIFIER;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Sets;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SnpEffAnnotationParserTest
{

    private VCFCodec codec;

    @Before
    public void setup()
    {
        codec = createTestCodec();
    }

    @Test
    public void testInitialSpliceBase()
    {
        assertEquals(4, SnpEffAnnotationParser.initialIndelSpliceBase(false, "c.5667+4_5667+9delAAGGGG"));
        assertEquals(3, SnpEffAnnotationParser.initialIndelSpliceBase(true, "c.1023+2_1023+3insA"));
        assertEquals(5, SnpEffAnnotationParser.initialIndelSpliceBase(false, "c.865+5dupT"));

        assertEquals(-1, SnpEffAnnotationParser.initialIndelSpliceBase(false, "c.2001+6_2001+40delAGCATTTATTTTAGGGGGAAAAAAAAAGTCGTACA"));
    }

    @Test
    public void testSnvSpliceDonorPlus5()
    {
        assertEffect("splice_region_variant",
                "G",
                "T",
                Strings.EMPTY,
                Strings.EMPTY,
                0,
                "splice_region_variant&intron_variant",
                "c.375+4G>T");
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "T", Strings.EMPTY, Strings.EMPTY, 0, "splice_region_variant&intron_variant", "c.375+5G>T");
        assertEffect(SPLICE_DONOR_CONSEQUENCE,
                "GG",
                "TT",
                Strings.EMPTY,
                Strings.EMPTY,
                0,
                "splice_region_variant&intron_variant",
                "c.375+5GG>TT");
        assertEffect("splice_region_variant",
                "G",
                "T",
                Strings.EMPTY,
                Strings.EMPTY,
                0,
                "splice_region_variant&intron_variant",
                "c.375+6G>T");

        assertEffect("intron_variant", "G", "T", Strings.EMPTY, Strings.EMPTY, 0, "intron_variant", "c.375+5G>T");
    }

    @Test
    public void testForwardStrandInsert()
    {
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TA", Strings.EMPTY, Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");

        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TA", "G", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TA", "GG", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect("splice_region_variant", "G", "TA", "GGG", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");

        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TA", Strings.EMPTY, "A", 2, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect("splice_region_variant", "G", "TA", Strings.EMPTY, "A", 3, "splice_region_variant", "c.1023+2_1023+3insA");

        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TA", Strings.EMPTY, "T", 24, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect("splice_region_variant", "G", "TTA", Strings.EMPTY, "TA", 24, "splice_region_variant", "c.1023+2_1023+3insTA");
    }

    @Test
    public void testReverseStrandInsert()
    {
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", Strings.EMPTY, Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");

        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", "G", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", "GG", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", "GGG", Strings.EMPTY, 0, "splice_region_variant", "c.1023+2_1023+3insA");

        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", Strings.EMPTY, "A", 2, "splice_region_variant", "c.1023+2_1023+3insA");
        assertEffect(SPLICE_DONOR_CONSEQUENCE, "G", "TT", Strings.EMPTY, "A", 3, "splice_region_variant", "c.1023+2_1023+3insA");
    }

    private void assertEffect(@NotNull String expectedStart, @NotNull String ref, @NotNull String alt, @NotNull String mh,
            @NotNull String repeat, int repeatCount, String effects, @NotNull String hgvs)
    {
        VariantContext context = create(ref, alt, mh, repeat, repeatCount, effects, hgvs);
        List<SnpEffAnnotation> annotations = SnpEffAnnotationParser.fromContext(context);
        assertEquals(1, annotations.size());
        assertTrue(annotations.get(0).effects().startsWith(expectedStart));
    }

    @NotNull
    private VariantContext create(@NotNull String ref, @NotNull String alt, @NotNull String mh, @NotNull String repeat, int repeatCount,
            @NotNull String effects, @NotNull String hgvs)
    {
        final String annotation = String.format("A|%s|LOW|GENE|GENEID|transcript|TRANSCIPTID|protein_coding|4/10|%s||||||", effects, hgvs);

        final StringJoiner infoJoiner = new StringJoiner(";").add(SNPEFF_IDENTIFIER + "=" + annotation);
        if(!mh.isEmpty())
        {
            infoJoiner.add(MICROHOMOLOGY_FLAG + "=" + mh);
        }
        if(repeatCount > 0)
        {
            infoJoiner.add(REPEAT_COUNT_FLAG + "=" + repeatCount);
            infoJoiner.add(REPEAT_SEQUENCE_FLAG + "=" + repeat);
        }

        final String line =
                String.format("15\t12345678\trs1;UCSC\t%s\t%s\t2\tPASS\t%s\tGT:AD:DP\t0/1:59,60:120", ref, alt, infoJoiner.toString());

        return codec.decode(line);
    }

    @NotNull
    private static VCFCodec createTestCodec()
    {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("SAMPLE"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
