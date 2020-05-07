package com.hartwig.hmftools.serve.vicc;

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

        LOGGER.info("Reading VICC json from {} with max entries 100", viccJsonPath);
        List<ViccEntry> viccEntries = ViccJsonReader.readViccKnowledgebaseJsonFileWithMaxEntries(viccJsonPath, 100);
        LOGGER.info("Read {} entries", viccEntries.size());

        String source = "oncokb";
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntriesSpecificKnowledgeBase = ViccJsonReader.readSingleKnowledgebase(viccJsonPath, source);
        LOGGER.info("Read {} entries", viccEntriesSpecificKnowledgeBase.size());
    }
}
