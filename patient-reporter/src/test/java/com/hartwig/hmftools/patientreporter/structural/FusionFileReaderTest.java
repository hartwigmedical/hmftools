package com.hartwig.hmftools.patientreporter.structural;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class FusionFileReaderTest {

    private static final String FUSION_FILE = Resources.getResource("test_run/svAnalysis/sample.linx.fusions.tsv").getPath();

    @Test
    public void canReadFromFile() throws IOException {
       assertEquals(1, FusionFileReader.fromFusionFile(FUSION_FILE).size());
    }
}
