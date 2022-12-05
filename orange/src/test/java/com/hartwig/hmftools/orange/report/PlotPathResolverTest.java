package com.hartwig.hmftools.orange.report;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class PlotPathResolverTest {

    @Test
    public void canResolvePlotPaths() {
        PlotPathResolver nullResolver = new PlotPathResolver(null);
        PlotPathResolver realResolver = new PlotPathResolver("/path/to/output");

        assertEquals("file", nullResolver.resolve("file"));
        assertEquals("/path/to/output/file", realResolver.resolve("file"));
    }
}