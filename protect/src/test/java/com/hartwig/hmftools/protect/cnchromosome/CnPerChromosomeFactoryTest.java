package com.hartwig.hmftools.protect.cnchromosome;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CnPerChromosomeFactoryTest {

    private static final String PURPLE_COPYNUMBER_TSV = Resources.getResource("cnchromosome/sample_purple.cnv.somatic.tsv").getPath();

    @Test
    public void extractCopyNumberPerChromosomeArm() throws IOException {
        CnPerChromosomeFactory.fromPurpleSomaticCopynumberTsv(PURPLE_COPYNUMBER_TSV);

    }

}