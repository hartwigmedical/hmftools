package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestUtils.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.INPUT_DIR;
import static com.hartwig.hmftools.chord.ChordTestUtils.TMP_OUTPUT_DIR;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class ChordRunnerTest
{
    @Before
    public void setup()
    {
        new File(TMP_OUTPUT_DIR).mkdir();
    }

    @After
    public void teardown() throws IOException
    {
        FileUtils.deleteDirectory(new File(TMP_OUTPUT_DIR));
    }

    @Test
    public void canRunChordOnEmptyVcfs()
    {
        String chordToolDir = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/chord/src/main/R/";

        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(List.of(EMPTY_SAMPLE))
                .purpleDir(INPUT_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .chordToolDir(chordToolDir)
                .build();

        ChordRunner runner = new ChordRunner(config);
        runner.run();
    }
}
