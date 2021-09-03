package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class AnnotatedVirusFileTest {

    private static final String SAMPLE_VIRUS_ANNOTATED_TSV = Resources.getResource("virus/sample.virus.annotated.tsv").getPath();

    @Test
    public void canReadSampleVirusAnnotatedTsv() throws IOException {
        List<AnnotatedVirus> annotatedVirusList = AnnotatedVirusFile.read(SAMPLE_VIRUS_ANNOTATED_TSV);
        assertEquals(1, annotatedVirusList.size());
    }

    @Test
    public void canConvertAnnotatedViruses() {
        AnnotatedVirus virus = VirusTestFactory.testAnnotatedVirusBuilder().build();

        List<String> lines = AnnotatedVirusFile.toLines(Lists.newArrayList(virus));
        List<AnnotatedVirus> converted = AnnotatedVirusFile.fromLines(lines);

        assertEquals(virus, converted.get(0));
    }
}