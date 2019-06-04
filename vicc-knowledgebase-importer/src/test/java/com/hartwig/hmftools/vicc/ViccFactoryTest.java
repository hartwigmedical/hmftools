package com.hartwig.hmftools.vicc;

import static org.junit.Assert.assertNotNull;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.junit.Ignore;
import org.junit.Test;

public class ViccFactoryTest {

    @Test
    @Ignore
    public void canReadJsonKnowledgeBaseFile() throws IOException {
        final String baseDir =
                System.getProperty("user.home") + File.separator + "hmf" + File.separator + "projects" + File.separator + "vicc";
        final String inputFile = baseDir + File.separator + "all.json";

        List<ViccEntry> viccEntries = ViccFactory.readViccKnowledgebaseJsonFile(inputFile);

        assertNotNull(viccEntries);
    }
}