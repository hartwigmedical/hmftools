package com.hartwig.hmftools.boggs.io;

import com.google.common.io.Resources;
import com.hartwig.hmftools.boggs.PatientData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagstatData2;
import com.hartwig.hmftools.boggs.flagstatreader.FlagstatParser2;
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

    private static class DummyFlagstatParser implements FlagstatParser2 {

        @NotNull
        public FlagstatData2 parse(@NotNull File file) throws IOException {
            // KODU (TODO): Mock this.
            return new FlagstatData2(file.getPath(), new FlagStatsBuilder().build(), new FlagStatsBuilder().build());
        }
    }
}
