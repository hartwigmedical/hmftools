package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Sets;

import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class SomaticVariantFactoryTest {

    private static final String SAMPLE = "sample";

    private SomaticVariantFactoryNew victim;
    private VCFCodec codec;

    private static final double EPSILON = 1.0e-10;

    @Before
    public void setup() {
        victim = new SomaticVariantFactoryNew();
        codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
    }

    @Test
    public void canReadSampleNameFromHeader() {
        final String header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";
        assertEquals("sample", SomaticVariantFactoryOld.sampleFromHeaderLine(header));
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
        final String missingAFLine = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9";
        final SomaticVariant missingAFVariant = SomaticVariantFactoryOld.fromVCFLine(missingAFLine);
        assertEquals(Double.NaN, missingAFVariant.alleleFrequency(), EPSILON);
        assertEquals(0, missingAFVariant.totalReadCount(), EPSILON);

        final String missingRefCovLine = "0 \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 0/1:60:113";
        final SomaticVariant missingRefCovVariant = SomaticVariantFactoryOld.fromVCFLine(missingRefCovLine);
        assertEquals(Double.NaN, missingRefCovVariant.alleleFrequency(), EPSILON);
        assertEquals(0, missingRefCovVariant.totalReadCount(), EPSILON);
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
        final SomaticVariant hasDBSNPOnly = SomaticVariantFactoryOld.fromVCFLine(dbsnpOnly);
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