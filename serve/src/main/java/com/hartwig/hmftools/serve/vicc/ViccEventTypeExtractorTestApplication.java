package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.vicc.eventtype.EventType;
import com.hartwig.hmftools.serve.vicc.eventtype.EventTypeAnalyzer;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;
import com.hartwig.hmftools.vicc.selection.ImmutableViccQuerySelection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccEventTypeExtractorTestApplication {

    private static final Logger LOGGER = LogManager.getLogger(ViccEventTypeExtractorTestApplication.class);

    public static void main(String[] args) throws IOException {
        String viccJsonPath = System.getProperty("user.home") + "/hmf/projects/vicc/all.json";

        ViccSource source = ViccSource.ONCOKB;
        LOGGER.info("Reading VICC json from {} with source '{}'", viccJsonPath, source);
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJsonPath,
                ImmutableViccQuerySelection.builder().addSourcesToFilterOn(source).build());
        LOGGER.info("Read {} entries", viccEntries.size());

        for (ViccEntry viccEntry : viccEntries) {
            LOGGER.info("Extracting event types for {}", viccEntry);

            List<EventType> eventTypes = EventTypeAnalyzer.determineEventTypes(viccEntry);
            LOGGER.info(" Extracted {} event types", eventTypes);
        }
    }
}
