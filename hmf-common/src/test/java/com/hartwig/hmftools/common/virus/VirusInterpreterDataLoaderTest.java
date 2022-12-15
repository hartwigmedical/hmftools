package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class VirusInterpreterDataLoaderTest {

    private static final String SAMPLE_VIRUS_ANNOTATED_TSV = Resources.getResource("virus/sample.virus.annotated.tsv").getPath();

    @Test
    public void canLoadVirusInterpreterData() throws IOException {
        VirusInterpreterData virusInterpreterData = VirusInterpreterDataLoader.load(SAMPLE_VIRUS_ANNOTATED_TSV);

        assertEquals(2, virusInterpreterData.allViruses().size());
        assertEquals(1, virusInterpreterData.reportableViruses().size());
    }
}