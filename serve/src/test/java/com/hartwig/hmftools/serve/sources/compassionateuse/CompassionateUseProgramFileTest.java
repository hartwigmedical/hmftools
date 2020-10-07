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

        CompassionateUseProgram compassionateUsePrograms1 = compassionateUsePrograms.get(0);
        assertEquals("AC", compassionateUsePrograms1.trialAcronymSite());
        assertEquals("-", compassionateUsePrograms1.variants());
        assertEquals("Source", compassionateUsePrograms1.source());
        assertEquals("Trametinib", compassionateUsePrograms1.drug());
        assertEquals("inhibitor", compassionateUsePrograms1.drugType());
        assertEquals("Lung", compassionateUsePrograms1.cancerType());
        assertEquals("B", compassionateUsePrograms1.level());
        assertEquals("Responsive", compassionateUsePrograms1.direction());
        assertEquals("link", compassionateUsePrograms1.link());
    }
}