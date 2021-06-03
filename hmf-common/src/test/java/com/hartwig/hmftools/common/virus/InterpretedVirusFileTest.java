package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class InterpretedVirusFileTest {

    private static final String SAMPLE_VIRUS_INTERPRETATION_TSV = Resources.getResource("virus/sample.virus_interpretation.tsv").getPath();

    @Test
    public void canReadSampleVirusInterpretationTsv() throws IOException {
        List<InterpretedVirus> interpretedVirusList = InterpretedVirusFile.read(SAMPLE_VIRUS_INTERPRETATION_TSV);
        assertEquals(2, interpretedVirusList.size());
    }

    @Test
    public void canConvertInterpretedViruses() {
        InterpretedVirus virus = VirusTestFactory.testInterpretedVirusBuilder().build();

        List<String> lines = InterpretedVirusFile.toLines(Lists.newArrayList(virus));
        List<InterpretedVirus> converted = InterpretedVirusFile.fromLines(lines);

        assertEquals(virus, converted.get(0));
    }
}