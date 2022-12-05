package com.hartwig.hmftools.orange.algo.plot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import org.junit.Test;

public class FileBasedPlotManagerTest {

    @Test
    public void canHandleNullFiles() throws IOException {
        FileBasedPlotManager plotManager = new FileBasedPlotManager("/path/to/output");
        assertNull(plotManager.processPlotFile(null));
    }

    @Test
    public void canExtractFileNames() {
        assertEquals("file.pdf", FileBasedPlotManager.extractFileName("/path/to/file.pdf"));
    }

    @Test
    public void canMakePathsRelative() {
        assertEquals("plot/picture.jpg", FileBasedPlotManager.relativePath("/full/path/to/plot/picture.jpg", "/full/path/to"));
        assertEquals("plot/picture.jpg", FileBasedPlotManager.relativePath("/full/path/to/plot/picture.jpg", "/full/path/to/"));
    }

    @Test (expected = IllegalStateException.class)
    public void crashOnInvalidRelativePathRequest() {
        FileBasedPlotManager.relativePath("/full/path/to/picture.jpg", "does not exist");
    }
}