package com.hartwig.hmftools.serve;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RefGenomeVersionTest {

    @Test
    public void canVersionFilePaths() {
        String path = "/this/is/my/path.vcf";
        assertEquals("/this/is/my/path.hg19.vcf", RefGenomeVersion.HG19.makeVersioned(path));

        String path2 = "file.testing.tsv";
        assertEquals("file.testing.hg19.tsv", RefGenomeVersion.HG19.makeVersioned(path2));
    }

    @Test(expected = IllegalStateException.class)
    public void cannotHandlePathsWithNoExtension() {
        RefGenomeVersion.HG19.makeVersioned("path");
    }
}