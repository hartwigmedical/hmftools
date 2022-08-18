package com.hartwig.hmftools.common.rose;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class RoseConclusionFileTest  {

    private static final String ROSE_TSV = Resources.getResource("rose/tumor_sample.rose.tsv").getPath();

    @Test
    public void canReadRoseSummaryFile() throws IOException {
        String roseSummary = RoseConclusionFile.read(ROSE_TSV);
        assertEquals("Melanoma sample <enter> - A <enter> - B <enter> C <enter> D <enter> ", roseSummary);
    }

}