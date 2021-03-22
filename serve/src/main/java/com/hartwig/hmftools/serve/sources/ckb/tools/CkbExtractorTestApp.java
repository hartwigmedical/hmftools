package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.JsonDatabaseToCkbEntryConverter;
import com.hartwig.hmftools.ckb.classification.CKBClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.CkbJsonReader;
import com.hartwig.hmftools.common.serve.classification.EventClassifier;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.EventClassifierFactory;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.ckb.CkBReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);
    private static final String FIELD_DELIMITER = "\t";

    public static void main(String[] args) throws IOException {
        String ckbDir =  "/data/common/dbs/ckb/210319_flex_dump";

        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        List<CkbEntry> ckbEntries = JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
      //  List<CkbEntry> curateCKBEntries = CkBReader.filterRelevantEntries(ckbEntries);


        EventClassifierConfig config = CKBClassificationConfig.build();
        EventClassifier classifier = EventClassifierFactory.buildClassifier(config);

        // TODO 1. Make sure every entry has correct event type.

        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("event").add("type").toString();
        lines.add(header);

        for (CkbEntry entry : ckbEntries) {
            String gene = entry.variants().get(0).gene().geneSymbol();
            String profileName = Strings.EMPTY;
            if (entry.profileName().split(" ")[0].equals(gene) && !entry.profileName().contains("-")) {
                profileName = entry.profileName().split(" ", 2)[1]; // remove gene name of profileName
            } else {
                profileName = entry.profileName();
            }

            EventType type = classifier.determineType(gene, profileName);
            if (entry.variants().size() > 1) {
                type = EventType.COMBINED;
            } else {
                LOGGER.info("Type of {} on {} is {}", profileName, gene, type);
                if (type == EventType.UNKNOWN) {
                    LOGGER.warn("Hier moet Lieke nog iets mee doen!!! {}", profileName);
                }

            }

            lines.add(new StringJoiner(FIELD_DELIMITER).add(profileName).add(type.toString()).toString());
        }
        Files.write(new File("/data/common/dbs/serve/pilot_output/events.tsv").toPath(), lines);


        // TODO 2. Make sure every event is extracted correctly using EventExtractor
//        EventExtractor extractor = EventExtractorFactory.create(config, ....);

        // TODO 3. Create ActionableEvents for all relevant entries.
    }
}
