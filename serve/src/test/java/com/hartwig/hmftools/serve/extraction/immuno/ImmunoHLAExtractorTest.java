package com.hartwig.hmftools.serve.extraction.immuno;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

public class ImmunoHLAExtractorTest {

    @Test
    public void canExtractImmunoHla() {
        ImmunoHLAExtractor immunoHLAExtractor = new ImmunoHLAExtractor();
        ImmunoHLA immunohLA = immunoHLAExtractor.extract(EventType.IMMUNO_HLA, "A*20");

        assertNotNull(immunohLA);
        assertEquals("A*20", immunohLA.immunoHLA());
    }
}