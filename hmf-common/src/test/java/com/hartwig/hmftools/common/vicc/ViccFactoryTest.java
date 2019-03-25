package com.hartwig.hmftools.common.vicc;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;
import org.junit.Ignore;

public class ViccFactoryTest {

    private static final String BRCA_FILE_PATH = Resources.getResource("vicc_knowledgebase").getPath();

    @Test
    @Ignore
    public void convertBRCA() throws IOException {
        ViccFactory.extractBRCAFile(BRCA_FILE_PATH + "/brca.json");
    }

    @Test
    @Ignore
    public void convertAll() throws IOException {
        ViccFactory.extractAllFile("all.json");
    }
}