package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class AnnotatedVirusFileV1Test {

    private static final String SAMPLE_VIRUS_ANNOTATED_TSV = Resources.getResource("virus/sample.virus.annotated.tsv").getPath();

    @Test
    public void canReadSampleVirusAnnotatedTsv() throws IOException {
        List<AnnotatedVirusV1> annotatedVirusList = AnnotatedVirusFileV1.read(SAMPLE_VIRUS_ANNOTATED_TSV);
        assertEquals(2, annotatedVirusList.size());
    }
}