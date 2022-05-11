package com.hartwig.hmftools.serve.extraction.immuno;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.junit.Test;

public class ImmunoHLAExtractorTest {

    @Test
    public void canExtractImmunoHLA() {
        ImmunoHLAExtractor immunoHLAExtractor = new ImmunoHLAExtractor();
        ImmunoHLA immunohLA = immunoHLAExtractor.extract(EventType.IMMUNO_HLA, "A*20");

        assertEquals("A*20", immunohLA.immunoHLA());
    }
}