package com.hartwig.hmftools.serve.sources.compassionateuse;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class CompassionateUseProgramFileTest {

    private static final String EXAMPLE_TSV = Resources.getResource("compassionateuse/example.tsv").getPath();

    @Test
    public void canReadTestCompassionateUseProgramFile() throws IOException {
        List<CompassionateUseProgram> compassionateUsePrograms = CompassionateUseProgramFile.read(EXAMPLE_TSV);

        assertEquals(1, compassionateUsePrograms.size());

        CompassionateUseProgram compassionateUseProgram = compassionateUsePrograms.get(0);
        assertEquals("AC", compassionateUseProgram.trialAcronymSite());
        assertEquals("-", compassionateUseProgram.variants());
        assertEquals("Source", compassionateUseProgram.source());
        assertEquals("Trametinib", compassionateUseProgram.drug());
        assertEquals("inhibitor", compassionateUseProgram.drugType());
        assertEquals("Lung", compassionateUseProgram.cancerType());
        assertEquals("B", compassionateUseProgram.level());
        assertEquals("Responsive", compassionateUseProgram.direction());
        assertEquals("link", compassionateUseProgram.link());
    }
}