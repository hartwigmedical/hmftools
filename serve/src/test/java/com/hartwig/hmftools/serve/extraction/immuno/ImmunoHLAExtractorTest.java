package com.hartwig.hmftools.serve.extraction.immuno;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.classification.EventType;

import org.junit.Test;

public class ImmunoHLAExtractorTest {

    @Test
    public void canExtractImmunoHLA() {
        ImmunoHLAExtractor immunoHLAExtractor = new ImmunoHLAExtractor();

        ImmunoHLA immunohLACorrect = immunoHLAExtractor.extract(EventType.IMMUNO_HLA, "A*20");
        assertEquals("A*20", immunohLACorrect.immunoHLA());

        ImmunoHLA immunohLALess = immunoHLAExtractor.extract(EventType.IMMUNO_HLA, "A*2");
        assertEquals("A*2", immunohLALess.immunoHLA());

        ImmunoHLA immunohLAMore = immunoHLAExtractor.extract(EventType.IMMUNO_HLA, "HLA-A*021:054");
        assertEquals("HLA-A*021:054", immunohLAMore.immunoHLA());
    }
}