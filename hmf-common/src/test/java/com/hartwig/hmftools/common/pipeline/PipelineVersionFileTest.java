package com.hartwig.hmftools.common.pipeline;

import static com.hartwig.hmftools.common.variant.VcfFileReader.parseToolVersion;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class PipelineVersionFileTest
{
    private static final String RESOURCE_DIR = Resources.getResource("pipeline").getPath();

    @Test
    public void canResolvePipelineVersionPre() throws IOException
    {
        String pipelineVersion = PipelineVersionFile.majorDotMinorVersion(RESOURCE_DIR + File.separator + "v5.7.pipeline.version");
        assertEquals("5.7", pipelineVersion);
    }

    @Test
    public void canResolvePipelineVersionPost() throws IOException
    {
        String pipelineVersion = PipelineVersionFile.majorDotMinorVersion(RESOURCE_DIR + File.separator + "v5.16.pipeline.version");
        assertEquals("5.16", pipelineVersion);
    }

    @Test
    public void canConvertToMajorDotMinor()
    {
        assertEquals("5.13", PipelineVersionFile.convertToMajorDotVersion("5.13.1234"));
        assertNull(PipelineVersionFile.convertToMajorDotVersion("local-snapshot"));
    }

    @Test
    public void testToolVcfVersion()
    {
        assertEquals(5.1, parseToolVersion("5.1"), 0.001);
        assertEquals(5.1, parseToolVersion("5.1.beta.1"), 0.001);
        assertEquals(5.1, parseToolVersion("5.1-abc"), 0.001);
        assertEquals(5, parseToolVersion("5.abc"), 0.001);
        assertNull(parseToolVersion("a5.1b"));
        assertNull(parseToolVersion(".1"));
    }
}