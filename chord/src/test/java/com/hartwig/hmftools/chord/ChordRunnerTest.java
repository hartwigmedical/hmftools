package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestUtils.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.HUMAN_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestUtils.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestUtils.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestUtils.TMP_OUTPUT_DIR;
import static com.hartwig.hmftools.chord.prep.ChordDataWriter.COHORT_FILE_PREFIX;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.chord.prep.ChordDataPrep;
import com.hartwig.hmftools.chord.prep.ChordDataWriter;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
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

    @Ignore
    @Test
    public void canPrepMultipleSamples()
    {
        Configurator.setRootLevel(Level.DEBUG);

        List<String> sampleIds = List.of(MINIMAL_SAMPLE, EMPTY_SAMPLE);

        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(sampleIds)
                .refGenomeFile(HUMAN_GENOME_FASTA)
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .threads(sampleIds.size())
                .build();

        ChordDataPrep prep = new ChordDataPrep(config);
        prep.run();

        File outputFile = new File(ChordDataWriter.formOutputFile(config, COHORT_FILE_PREFIX));
        assertTrue(outputFile.exists());
    }

    @Test
    public void canRunChordOnEmptyVcfs()
    {
        String chordToolDir = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/chord/src/main/R/";

        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(List.of(EMPTY_SAMPLE))
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .chordToolDir(chordToolDir)
                .build();

        ChordRunner runner = new ChordRunner(config);
        runner.run();
    }
}
