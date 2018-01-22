package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SomaticVariantFactoryTest {

    private static final String SAMPLE = "sample";
    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";

    private SomaticVariantFactory victim;
    private VCFCodec codec;

    private static final double EPSILON = 1.0e-10;

    @Before
    public void setup() {
        victim = new SomaticVariantFactory();
        codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
    }

    @Test
    public void canLoadSomaticVCFFromBasePath() throws IOException, HartwigException {
        final List<SomaticVariant> variants =  new SomaticVariantFactory().fromVCFFile("sample", VARIANT_PATH, SOMATIC_EXTENSION);
        assertEquals(3, variants.size());
    }

    @Test
    public void canLoadSomaticVCFFromFile() throws IOException, HartwigException {
        final String file = VARIANT_PATH + File.separator + SOMATIC_EXTENSION;
        final List<SomaticVariant> variants =  new SomaticVariantFactory().fromVCFFile("sample", file);
        assertEquals(3, variants.size());
    }

    @Test
    public void ignoreZeroReadCount() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tset=varscan-freebayes;\tGT:AD:DP\t0/1:0,0:121";
        assertFalse(victim.createVariant(SAMPLE, codec.decode(line)).isPresent());
    }

    @Test
    public void canReadAndWriteCorrectSomaticVariant() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tset=varscan-freebayes;\tGT:AD:DP\t0/1:60,60:121";

        final SomaticVariant variant = victim.createVariant(SAMPLE, codec.decode(line)).get();
        assertEquals("15", variant.chromosome());
        assertEquals(12345678, variant.position());
        assertEquals(VariantType.SNP, variant.type());
        assertTrue(variant.isDBSNP());
        assertFalse(variant.isCOSMIC());

        assertEquals("C", variant.ref());
        assertEquals("A,G", variant.alt());

        assertEquals(2, variant.callerCount());
        assertTrue(variant.callers().contains(SomaticVariantConstants.FREEBAYES));
        assertTrue(variant.callers().contains(SomaticVariantConstants.VARSCAN));
        assertFalse(variant.callers().contains(SomaticVariantConstants.STRELKA));
        assertFalse(variant.callers().contains(SomaticVariantConstants.MUTECT));

        assertEquals(0.5, variant.alleleFrequency(), EPSILON);
        assertEquals(120, variant.totalReadCount(), EPSILON);
    }

    @Test
    public void incorrectSampleFieldYieldsMissingReadCounts() {
        final String missingAFLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tset=varscan-freebayes;\tGT:DP\t0/1:21";

        final Optional<SomaticVariant> missingAFVariant = victim.createVariant(SAMPLE, codec.decode(missingAFLine));
        assertFalse(missingAFVariant.isPresent());

        final String missingRefCovLine = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tset=varscan-freebayes;\tGT:AD:DP\t0/1:60:121";
        final Optional<SomaticVariant> missingRefCovVariant = victim.createVariant(SAMPLE, codec.decode(missingRefCovLine));
        assertFalse(missingRefCovVariant.isPresent());
    }

    @Test
    public void recognizeFilterInVarscan() {
        final String line = "15\t12345678\trs1;UCSC\tC\tA,G\t2\tPASS\tset=freebayes-filterInVarscan;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant variant = victim.createVariant(SAMPLE, codec.decode(line)).get();

        assertEquals(1, variant.callerCount());
        assertFalse(variant.callers().contains(SomaticVariantConstants.VARSCAN));
    }

    @Test
    public void correctDBSNPAndCOSMIC() {
        final String both = "15\t12345678\trs1;COSM2\tC\tA,G\t2\tPASS\tset=freebayes-filterInVarscan;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasBoth = victim.createVariant(SAMPLE, codec.decode(both)).get();
        assertTrue(hasBoth.isDBSNP());
        assertTrue(hasBoth.isCOSMIC());
        assertEquals("COSM2", hasBoth.cosmicID());

        final String dbsnpOnly = "15\t12345678\trs1\tC\tA,G\t2\tPASS\tset=freebayes-filterInVarscan;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasDBSNPOnly = victim.createVariant(SAMPLE, codec.decode(dbsnpOnly)).get();
        assertTrue(hasDBSNPOnly.isDBSNP());
        assertFalse(hasDBSNPOnly.isCOSMIC());

        final String cosmicOnly = "15\t12345678\tCOSM2\tC\tA,G\t2\tPASS\tset=freebayes-filterInVarscan;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasCOSMICOnly = victim.createVariant(SAMPLE, codec.decode(cosmicOnly)).get();
        assertFalse(hasCOSMICOnly.isDBSNP());
        assertTrue(hasCOSMICOnly.isCOSMIC());
        assertEquals("COSM2", hasCOSMICOnly.cosmicID());

        final String none = "15\t12345678\tID\tC\tA,G\t2\tPASS\tset=freebayes-filterInVarscan;\tGT:AD:DP\t0/1:60,60:121";
        final SomaticVariant hasNone = victim.createVariant(SAMPLE, codec.decode(none)).get();
        assertFalse(hasNone.isDBSNP());
        assertFalse(hasNone.isCOSMIC());
    }

}