package com.hartwig.hmftools.common.genome.refgenome;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RefGenomeVersionTest {

    @Test
    public void canVersionChromosomes() {
        String chr37 = "10";
        String chr19 = "chr10";
        String chr38 = "chr10";

        assertEquals(chr37, V37.versionedChromosome(chr19));
        assertEquals(chr37, V37.versionedChromosome(chr37));
        assertEquals(chr37, V37.versionedChromosome(chr38));

        assertEquals(chr38, V38.versionedChromosome(chr19));
        assertEquals(chr38, V38.versionedChromosome(chr37));
        assertEquals(chr38, V38.versionedChromosome(chr38));
    }

    @Test
    public void canVersionFilePaths() {
        String path = "/this/is/my/path.vcf";
        assertEquals("/this/is/my/path.37.vcf", RefGenomeVersion.V37.addVersionToFilePath(path));

        String path2 = "file.testing.tsv";
        assertEquals("file.testing.37.tsv", RefGenomeVersion.V37.addVersionToFilePath(path2));

        String path3 = "file.vcf.gz";
        assertEquals("file.37.vcf.gz", RefGenomeVersion.V37.addVersionToFilePath(path3));
    }

    @Test(expected = IllegalStateException.class)
    public void cannotHandlePathsWithNoExtension() {
        RefGenomeVersion.V37.addVersionToFilePath("path");
    }

    @Test(expected = IllegalStateException.class)
    public void cannotHandlePathWithJustGzipExtension() {
        RefGenomeVersion.V37.addVersionToFilePath("path.gz");
    }
}