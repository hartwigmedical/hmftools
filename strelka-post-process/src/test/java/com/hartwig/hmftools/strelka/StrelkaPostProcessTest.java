package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.StrelkaPostProcess.allelicFrequency;
import static com.hartwig.hmftools.strelka.StrelkaPostProcess.qualityScore;
import static com.hartwig.hmftools.strelka.StrelkaPostProcess.variantGenomePosition;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;

import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class StrelkaPostProcessTest {
    private static final File VCF_FILE = new File(Resources.getResource("variants.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final String HIGH_CONF_BED_PATH = Resources.getResource("high_conf.bed").getPath();
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    private Slicer hcBed;
    private StrelkaPostProcess victim;

    @Before
    public void setup() throws IOException {
        hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        victim = new StrelkaPostProcess(hcBed);
    }

    @Test
    public void canSimplifySNP() {
        final VariantContext snpVariant = VARIANTS.get(0);
        assertTrue(snpVariant.isSNP());
        assertArrayEquals(new int[] { 39, 7 }, StrelkaPostProcess.getAD(snpVariant));
        assertEquals(53, StrelkaPostProcess.getDP(snpVariant));
    }

    @Test
    public void canSimplifyIndel() {
        final VariantContext indelVariant = VARIANTS.get(1);
        assertTrue(indelVariant.isIndel());
        assertArrayEquals(new int[] { 57, 9 }, StrelkaPostProcess.getAD(indelVariant));
        assertEquals(66, StrelkaPostProcess.getDP(indelVariant));
    }

    @Test
    public void canSimplifyMultipleAltVariant() {
        final VariantContext snpVariant = VARIANTS.get(2);
        assertTrue(snpVariant.isSNP());
        assertArrayEquals(new int[] { 0, 4, 1 }, StrelkaPostProcess.getAD(snpVariant));
        assertEquals(5, StrelkaPostProcess.getDP(snpVariant));
    }

    @Test
    public void canSimplifyMultipleAltIndelVariant() {
        final VariantContext indelVariant = VARIANTS.get(15);
        assertTrue(indelVariant.isIndel());
        assertArrayEquals(new int[] { 87, 3, 3 }, StrelkaPostProcess.getAD(indelVariant));
        assertEquals(100, StrelkaPostProcess.getDP(indelVariant));
    }

    @Test
    public void canSimplifyVariantWithoutRead() {
        final VariantContext snpVariant = VARIANTS.get(3);
        assertTrue(snpVariant.isSNP());
        assertArrayEquals(new int[] { 39, 0 }, StrelkaPostProcess.getAD(snpVariant));
        assertEquals(53, StrelkaPostProcess.getDP(snpVariant));
    }

    @Test
    public void canSimplifyVariantWithDotRead() {
        final VariantContext snpVariant = VARIANTS.get(4);
        assertTrue(snpVariant.isSNP());
        assertArrayEquals(new int[] { 39, 0 }, StrelkaPostProcess.getAD(snpVariant));
        assertEquals(53, StrelkaPostProcess.getDP(snpVariant));
    }

    @Test
    public void doesNotFilterMultipleAltVariant() {
        final VariantContext snpVariant = VARIANTS.get(2);
        assertTrue(victim.test(snpVariant));
    }

    @Test
    public void filtersVariantWithNoCall() {
        final VariantContext snpVariant = VARIANTS.get(5);
        assertFalse(victim.test(snpVariant));
    }

    @Test
    public void filtersHCRegionSNPBelowThreshold() {
        final VariantContext snpVariant = VARIANTS.get(6);
        assertEquals(0.5, allelicFrequency(snpVariant), 0.001);
        assertEquals(2, qualityScore(snpVariant));
        assertTrue(hcBed.includes(variantGenomePosition(snpVariant)));
        assertFalse(victim.test(snpVariant));
    }

    @Test
    public void doesNotFilterHCRegionSNPAboveThreshold() {
        final VariantContext snpVariant = VARIANTS.get(7);
        assertEquals(0.5, allelicFrequency(snpVariant), 0.001);
        assertEquals(3, qualityScore(snpVariant));
        assertTrue(hcBed.includes(variantGenomePosition(snpVariant)));
        assertTrue(victim.test(snpVariant));
    }

    @Test
    public void filtersLCRegionSNPBelowThreshold() {
        final VariantContext snpVariant = VARIANTS.get(8);
        assertEquals(0.1, allelicFrequency(snpVariant), 0.001);
        assertEquals(20, qualityScore(snpVariant));
        assertTrue(!hcBed.includes(variantGenomePosition(snpVariant)));
        assertFalse(victim.test(snpVariant));
    }

    @Test
    public void doesNotFilterLCRegionSNPAboveThreshold() {
        final VariantContext snpVariant = VARIANTS.get(9);
        assertEquals(0.11, allelicFrequency(snpVariant), 0.01);
        assertEquals(21, qualityScore(snpVariant));
        assertTrue(!hcBed.includes(variantGenomePosition(snpVariant)));
        assertTrue(victim.test(snpVariant));
    }

    @Test
    public void filtersHCRegionIndelBelowThreshold() {
        final VariantContext indelVariant = VARIANTS.get(10);
        assertEquals(0.5, allelicFrequency(indelVariant), 0.001);
        assertEquals(2, qualityScore(indelVariant));
        assertTrue(hcBed.includes(variantGenomePosition(indelVariant)));
        assertFalse(victim.test(indelVariant));
    }

    @Test
    public void doesNotFilterHCRegionIndelBelowThreshold() {
        final VariantContext indelVariant = VARIANTS.get(11);
        assertEquals(0.5, allelicFrequency(indelVariant), 0.001);
        assertEquals(3, qualityScore(indelVariant));
        assertTrue(hcBed.includes(variantGenomePosition(indelVariant)));
        assertTrue(victim.test(indelVariant));
    }

    @Test
    public void filtersLCRegionIndelBelowThreshold() {
        final VariantContext indelVariant = VARIANTS.get(12);
        assertEquals(0.1, allelicFrequency(indelVariant), 0.001);
        assertEquals(20, qualityScore(indelVariant));
        assertTrue(!hcBed.includes(variantGenomePosition(indelVariant)));
        assertFalse(victim.test(indelVariant));
    }

    @Test
    public void doesNotFilterLCRegionIndelAboveThreshold() {
        final VariantContext indelVariant = VARIANTS.get(13);
        assertEquals(0.11, allelicFrequency(indelVariant), 0.01);
        assertEquals(21, qualityScore(indelVariant));
        assertTrue(!hcBed.includes(variantGenomePosition(indelVariant)));
        assertTrue(victim.test(indelVariant));
    }

    @Test
    public void doesNotFilterHotspot() {
        final VariantContext indelVariant = VARIANTS.get(14);
        assertEquals(0.1, allelicFrequency(indelVariant), 0.001);
        assertEquals(20, qualityScore(indelVariant));
        assertTrue(indelVariant.hasAttribute("HOTSPOT"));
        assertTrue(!hcBed.includes(variantGenomePosition(indelVariant)));
        assertTrue(victim.test(indelVariant));
    }
}
