package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordTestDataPaths.EMPTY_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.DUMMY_GENOME_FASTA;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.INPUT_VCF_DIR;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.MINIMAL_SAMPLE;
import static com.hartwig.hmftools.chord.ChordTestDataPaths.TMP_OUTPUT_DIR;

import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import org.junit.Ignore;

@Ignore
public class ChordApplicationTest
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
        public void canPrepAndPredictOneSample()
        {
            Configurator.setRootLevel(Level.DEBUG);

            ChordConfig config = new ChordConfig.Builder()
                    .sampleIds(MINIMAL_SAMPLE)
                    .refGenomeFile(DUMMY_GENOME_FASTA)
                    .purpleDir(INPUT_VCF_DIR)
                    .outputDir(TMP_OUTPUT_DIR)
                    .build();

            ChordApplication chordApplication = new ChordApplication(config);
            chordApplication.run();

            File mutationContextsFile = new File(ChordOutput.mutationContextsFile(config));
            File predictionsFile = new File(ChordOutput.predictionsFile(config));

            assertTrue(mutationContextsFile.exists());
            assertTrue(predictionsFile.exists());
        }

    @Test
    public void canPrepAndPredictMultipleSamples()
    {
        Configurator.setRootLevel(Level.DEBUG);

        List<String> sampleIds = List.of(MINIMAL_SAMPLE, EMPTY_SAMPLE);

        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(sampleIds)
                .refGenomeFile(DUMMY_GENOME_FASTA)
                .purpleDir(INPUT_VCF_DIR)
                .outputDir(TMP_OUTPUT_DIR)
                .threads(sampleIds.size())
                .build();

        ChordApplication chordApplication = new ChordApplication(config);
        chordApplication.run();

        File mutationContextsFile = new File(ChordOutput.mutationContextsFile(config));
        File predictionsFile = new File(ChordOutput.predictionsFile(config));

        assertTrue(mutationContextsFile.exists());
        assertTrue(predictionsFile.exists());
    }
}
