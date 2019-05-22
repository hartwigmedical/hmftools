package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.IOException;

import org.junit.Ignore;
import org.junit.Test;

public class ViccFactoryTest {

    @Test
    @Ignore
    public void convertAll() throws IOException {
        final String baseDir = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
        final String inputFile = baseDir + File.separator + "all.json";
        final String outputCsvFileName = baseDir + File.separator + "all.csv";

        ViccFactory.extractAllFileSpecificFields(inputFile, outputCsvFileName);
    }
}