package com.hartwig.hmftools.sullivan;

import com.google.common.io.Resources;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;
import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SullivanAlgoTest {

    @Test
    public void sullivanAlgoRunsOnRealTestData() {
        URL originalFastqURL = Resources.getResource("fastq/original.fastq");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated.fastq");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        assertTrue(SullivanAlgo.runSullivanAlgo(originalFastqPath, recreatedFastqPath, null, false, 100000));
    }

    @Test
    public void sullivanAlgoRunsOnPartialData() {
        URL originalFastqURL = Resources.getResource("fastq/original.fastq");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated.fastq");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        // KODU: Assume that test files have 25k records.
        assertTrue(SullivanAlgo.runSullivanAlgo(originalFastqPath, recreatedFastqPath, null, false, 1000));
    }

    @Test
    public void convertsOriginalToRecreatedFileNames() {
        String original = "Set1GIAB12878_AHJKJVCCXX_S1_L001_R1_001.fastq.gz";
        String recreated = "Set1GIAB12878_AHJKJVCCXX_S1_L001_001_1.fastq";

        assertEquals(recreated, SullivanAlgo.fromOriginalToRecreatedFileName(original));
    }
}
