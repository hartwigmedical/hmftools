package com.hartwig.hmftools.patientdb.clinical.readers;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientdb.clinical.context.RunContext;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RunsJsonReaderTest {

    @SuppressWarnings("UnstableApiUsage")
    private static final String JSON_DIR = Resources.getResource("context" + File.separator + "RunJson").getPath();

    @Test(expected = IllegalArgumentException.class)
    public void throwsIllegalArgumentIfFileDoesntExist() {
        RunsJsonReader.extractRunContexts(new File("/does/not/exist"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void throwsIllegalArgumentOnEmptyFile() {
        RunsJsonReader.extractRunContexts(testJson("empty.json"));
    }

    @Test
    public void readsRunContextsFromJsonFile() {
        List<RunContext> runContexts = RunsJsonReader.extractRunContexts(testJson("two_somatic.json"));
        assertEquals(2, runContexts.size());
    }

    @NotNull
    private File testJson(final String filename) {
        return new File(JSON_DIR + File.separator + filename);
    }
}