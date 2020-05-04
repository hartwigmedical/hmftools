package com.hartwig.hmftools.knowledgebasegenerator.vicc;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.knowledgebasegenerator.vicc.eventtype.EventType;
import com.hartwig.hmftools.knowledgebasegenerator.vicc.eventtype.EventTypeAnalyzer;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccEventTypeExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccEventTypeExtractorTestApplication.class);

    public static void main(String[] args) throws IOException {
        String viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

        String source = "oncokb";
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntries = ViccJsonReader.readSingleKnowledgebase(viccJsonPath, source);
        LOGGER.info("Read {} entries", viccEntries.size());

        for (ViccEntry viccEntry : viccEntries) {
            LOGGER.info("Extracting event types for {}", viccEntry);

            List<EventType> eventTypes = EventTypeAnalyzer.determineEventTypes(viccEntry);
            LOGGER.info(" Extracted {} event types", eventTypes);

        }
    }
}
