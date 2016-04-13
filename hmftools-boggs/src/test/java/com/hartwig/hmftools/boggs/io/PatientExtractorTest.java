package com.hartwig.hmftools.boggs.io;

import com.google.common.io.Resources;
import com.hartwig.hmftools.boggs.PatientData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatParser;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatsBuilder;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import static org.junit.Assert.assertNotNull;

public class PatientExtractorTest {

    @Test
    public void canProcessRunDirectoryStructure() throws IOException {
        URL runDirURL = Resources.getResource("rundir");

        PatientExtractor extractor = new PatientExtractor(new DummyFlagstatParser());
        PatientData patient = extractor.extractFromRunDirectory(runDirURL.getPath());

        assertNotNull(patient);
    }

    private static class DummyFlagstatParser implements FlagStatParser {

        @NotNull
        public FlagStatData parse(@NotNull File file) throws IOException {
            // KODU (TODO): Mock this.
            return new FlagStatData(file.getPath(), new FlagStatsBuilder().build(), new FlagStatsBuilder().build());
        }
    }
}
