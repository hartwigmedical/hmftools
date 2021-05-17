package com.hartwig.hmftools.protect.cnchromosome;

import static org.junit.Assert.*;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CnPerChromosomeFactoryTest {

    @Test
    public void extarctCopyNumberperChromosomeArm() throws IOException {
        CnPerChromosomeFactory.fromPurpleSomaticCopynumberTsv(System.getProperty("user.home") + "/hmf/tmp/test.purple.cnv.somatic.tsv");

    }

}