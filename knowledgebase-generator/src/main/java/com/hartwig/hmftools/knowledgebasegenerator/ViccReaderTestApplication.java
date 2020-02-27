package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccReaderTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccReaderTestApplication.class);

    public static void main(String[] args) throws IOException {
        String viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

        LOGGER.info("Reading VICC json from {}", viccJsonPath);
       // List<ViccEntry> viccEntries = ViccJsonReader.readViccKnowledgebaseJsonFileWithMaxEntries(viccJsonPath, 100);
        List<ViccEntry> viccEntriesSpecificKnowledgeBase = ViccJsonReader.readViccKnowledgebaseJsonFileWithSpecificKnowledgebase(viccJsonPath, "jax_trials");

        LOGGER.info("Read {} entries", viccEntriesSpecificKnowledgeBase.size());
    }
}
