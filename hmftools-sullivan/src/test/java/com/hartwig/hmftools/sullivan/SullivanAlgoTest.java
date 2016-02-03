package com.hartwig.hmftools.sullivan;

import com.google.common.io.Resources;
import org.junit.Test;

import java.net.URL;

public class SullivanAlgoTest {

    @Test
    public void sullivanAlgoRunsOnRealTestData() {
        URL originalFastqURL = Resources.getResource("fastq/original.fastq");
        String originalFastqPath = originalFastqURL.getPath();

        URL recreatedFastqURL = Resources.getResource("fastq/recreated.fastq");
        String recreatedFastqPath = recreatedFastqURL.getPath();

        SullivanAlgo.runSullivanAlgo(originalFastqPath, recreatedFastqPath);
    }
}
