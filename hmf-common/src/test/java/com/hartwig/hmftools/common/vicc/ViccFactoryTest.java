package com.hartwig.hmftools.common.vicc;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;

public class ViccFactoryTest {

    private static final String BRCA_FILE_PATH = Resources.getResource("vicc_knowledgebase").getPath();

    @Test
    @Ignore
    public void convertBRCA() throws IOException {
        ViccFactory.extractBRCAFile(BRCA_FILE_PATH + "/brca.json");
    }

}