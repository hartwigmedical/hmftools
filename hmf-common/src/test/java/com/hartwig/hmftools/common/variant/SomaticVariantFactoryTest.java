package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SomaticVariantFactoryTest
{
    private static final String SAMPLE = "sample";
    private static final String SOMATIC_VARIANT_FILE = Resources.getResource("variant/somatics.vcf").getPath();

    private static final double EPSILON = 1.0e-10;

    private SomaticVariantFactory victim;
    private VCFCodec codec;

    @Before
    public void setup()
    {
        victim = new SomaticVariantFactory();
        codec = createTestCodec();
    }

    @NotNull
    private static VCFCodec createTestCodec()
    {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @Test
    public void canLoadSomaticSimpleVCFFromFile() throws IOException
    {
        final List<SomaticVariant> unfiltered = new SomaticVariantFactory().fromVCFFile("sample", SOMATIC_VARIANT_FILE);
        assertEquals(3, unfiltered.size());

        final List<SomaticVariant> filtered = SomaticVariantFactory.passOnlyInstance().fromVCFFile("sample", SOMATIC_VARIANT_FILE);
        assertEquals(2, filtered.size());
    }

    @Test
    public void canReadCorrectSomaticVariant()
    {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:59,60:120";

        final SomaticVariant variant = assertedGet(victim.createVariant(SAMPLE, codec.decode(line)));
        assertEquals("15", variant.chromosome());
        assertEquals(12345678, variant.position());
        assertEquals(VariantType.SNP, variant.type());

        assertEquals("C", variant.ref());
        assertEquals("A,G", variant.alt());

        assertEquals(0.5, variant.allelicDepth().alleleFrequency(), EPSILON);
        assertEquals(120, variant.allelicDepth().TotalReadCount, EPSILON);
    }

    @Test
    public void handleZeroReadCount()
    {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:0,0:0";
        final Optional<SomaticVariant> missingAFVariant = victim.createVariant(SAMPLE, codec.decode(line));
        assertFalse(missingAFVariant.isPresent());
    }

    @Test
    public void incorrectSampleFieldYieldsMissingReadCounts()
    {
        final String missingAFLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:DP\t0/1:21";

        final Optional<SomaticVariant> missingAFVariant = victim.createVariant(SAMPLE, codec.decode(missingAFLine));
        assertFalse(missingAFVariant.isPresent());

        final String missingRefCovLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tinfo;\tGT:AD:DP\t0/1:60:121";
        final Optional<SomaticVariant> missingRefCovVariant = victim.createVariant(SAMPLE, codec.decode(missingRefCovLine));
        assertFalse(missingRefCovVariant.isPresent());
    }

    @NotNull
    private static SomaticVariant assertedGet(@NotNull Optional<SomaticVariant> variant)
    {
        assert variant.isPresent();
        return variant.get();
    }
}