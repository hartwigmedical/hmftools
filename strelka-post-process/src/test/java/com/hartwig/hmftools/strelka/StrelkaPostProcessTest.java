package com.hartwig.hmftools.strelka;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariant;
import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariantFactory;

import org.junit.Test;

public class StrelkaPostProcessTest {
    private static final File VCF_FILE = new File(Resources.getResource("variants.vcf").getPath());
    private static final String HIGH_CONF_BED_PATH = Resources.getResource("high_conf.bed").getPath();

    @Test
    public void canSimplifySNP() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(0));
        assertEquals(VariantType.SNP, snpVariant.type());
        assertEquals("0/1:39,7:53", StrelkaPostProcessApplication.simplifyTumorDataString(snpVariant));
    }

    @Test
    public void canSimplifyIndel() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final StrelkaSomaticVariant indelVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(1));
        assertEquals(VariantType.INDEL, indelVariant.type());
        assertEquals("0/1:57,9:66", StrelkaPostProcessApplication.simplifyTumorDataString(indelVariant));
    }

    @Test
    public void canSimplifyMultipleAltVariant() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(2));
        assertEquals(VariantType.SNP, snpVariant.type());
        assertEquals("0/1:0,4,1:5", StrelkaPostProcessApplication.simplifyTumorDataString(snpVariant));
    }

    @Test
    public void canSimplifyVariantWithoutRead() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(3));
        assertEquals(VariantType.SNP, snpVariant.type());
        assertEquals("0/1:39,0:53", StrelkaPostProcessApplication.simplifyTumorDataString(snpVariant));
    }

    @Test
    public void canSimplifyVariantWithDotRead() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(4));
        assertEquals(VariantType.SNP, snpVariant.type());
        assertEquals("0/1:39,0:53", StrelkaPostProcessApplication.simplifyTumorDataString(snpVariant));
    }

    @Test
    public void doesNotFilterMultipleAltVariant() throws IOException, EmptyFileException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(2));
        assertTrue(StrelkaPostProcessApplication.checkVariant(snpVariant, hcBed));
    }

    @Test
    public void filtersHCRegionSNPBelowThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(6));
        assertEquals(0.5, snpVariant.allelicFrequency(), 0.001);
        assertEquals(2, snpVariant.qualityScore());
        assertTrue(hcBed.includes(snpVariant));
        assertFalse(StrelkaPostProcessApplication.checkVariant(snpVariant, hcBed));
    }

    @Test
    public void doesNotFilterHCRegionSNPAboveThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(7));
        assertEquals(0.5, snpVariant.allelicFrequency(), 0.001);
        assertEquals(3, snpVariant.qualityScore());
        assertTrue(hcBed.includes(snpVariant));
        assertTrue(StrelkaPostProcessApplication.checkVariant(snpVariant, hcBed));
    }

    @Test
    public void filtersLCRegionSNPBelowThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(8));
        assertEquals(0.1, snpVariant.allelicFrequency(), 0.001);
        assertEquals(20, snpVariant.qualityScore());
        assertTrue(!hcBed.includes(snpVariant));
        assertFalse(StrelkaPostProcessApplication.checkVariant(snpVariant, hcBed));
    }

    @Test
    public void doesNotFilterLCRegionSNPAboveThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant snpVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(9));
        assertEquals(0.11, snpVariant.allelicFrequency(), 0.01);
        assertEquals(21, snpVariant.qualityScore());
        assertTrue(!hcBed.includes(snpVariant));
        assertTrue(StrelkaPostProcessApplication.checkVariant(snpVariant, hcBed));
    }

    @Test
    public void filtersHCRegionIndelBelowThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant indelVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(10));
        assertEquals(0.5, indelVariant.allelicFrequency(), 0.001);
        assertEquals(2, indelVariant.qualityScore());
        assertTrue(hcBed.includes(indelVariant));
        assertFalse(StrelkaPostProcessApplication.checkVariant(indelVariant, hcBed));
    }

    @Test
    public void doesNotFilterHCRegionIndelBelowThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant indelVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(11));
        assertEquals(0.5, indelVariant.allelicFrequency(), 0.001);
        assertEquals(3, indelVariant.qualityScore());
        assertTrue(hcBed.includes(indelVariant));
        assertTrue(StrelkaPostProcessApplication.checkVariant(indelVariant, hcBed));
    }

    @Test
    public void filtersLCRegionIndelBelowThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant indelVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(12));
        assertEquals(0.1, indelVariant.allelicFrequency(), 0.001);
        assertEquals(20, indelVariant.qualityScore());
        assertTrue(!hcBed.includes(indelVariant));
        assertFalse(StrelkaPostProcessApplication.checkVariant(indelVariant, hcBed));
    }

    @Test
    public void doesNotFilterLCRegionIndelAboveThreshold() throws IOException, HartwigException {
        final List<String> vcfLines = FileReader.build().readLines(VCF_FILE.toPath());
        final Slicer hcBed = SlicerFactory.fromBedFile(HIGH_CONF_BED_PATH);
        final StrelkaSomaticVariant indelVariant = StrelkaSomaticVariantFactory.fromVCFLine(vcfLines.get(13));
        assertEquals(0.11, indelVariant.allelicFrequency(), 0.01);
        assertEquals(21, indelVariant.qualityScore());
        assertTrue(!hcBed.includes(indelVariant));
        assertTrue(StrelkaPostProcessApplication.checkVariant(indelVariant, hcBed));
    }
}
