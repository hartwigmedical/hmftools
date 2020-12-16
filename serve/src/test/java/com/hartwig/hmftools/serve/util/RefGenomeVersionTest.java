package com.hartwig.hmftools.serve.util;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RefGenomeVersionTest {

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