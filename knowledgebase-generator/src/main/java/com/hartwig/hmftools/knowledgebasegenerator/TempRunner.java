package com.hartwig.hmftools.knowledgebasegenerator;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TempRunner {

    private static final Logger LOGGER = LogManager.getLogger(TempRunner.class);

    public static void main(String[] args) throws IOException {
        // Need a decent memory to run this. -Xmx8G does the trick under VM options.
        String viccJson = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

        LOGGER.info("Loading VICC entries from {}", viccJson);
        List<ViccEntry> viccEntries = ViccJsonReader.readViccKnowledgebaseJsonFile(viccJson);
        LOGGER.info("Read {} VICC entries from {}", viccEntries.size(), viccJson);
    }
}
