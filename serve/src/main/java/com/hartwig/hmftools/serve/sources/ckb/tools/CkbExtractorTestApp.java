package com.hartwig.hmftools.serve.sources.ckb.tools;

import java.io.IOException;
import java.util.List;

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

public class CkbExtractorTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractorTestApp.class);

    public static void main(String[] args) throws IOException {
        String ckbDir =  "/data/common/dbs/ckb/210312_flex_dump";

        CkbJsonDatabase ckbJsonDatabase = CkbJsonReader.read(ckbDir);
        List<CkbEntry> ckbEntries = JsonDatabaseToCkbEntryConverter.convert(ckbJsonDatabase);
      //  List<CkbEntry> curateCKBEntries = CkBReader.filterRelevantEntries(ckbEntries);


        EventClassifierConfig config = CKBClassificationConfig.build();
        EventClassifier classifier = EventClassifierFactory.buildClassifier(config);

        // TODO 1. Make sure every entry has correct event type.
        for (CkbEntry entry : ckbEntries) {
            String gene = entry.variants().get(0).gene().geneSymbol();
            EventType type = classifier.determineType(gene, entry.profileName());
            LOGGER.info("Type of {} on {} is {}", entry.profileName(), gene, type);
            if (type == EventType.UNKNOWN) {
                LOGGER.warn("Hier moet Lieke nog iets mee doen!!!");
            }
        }

        // TODO 2. Make sure every event is extracted correctly using EventExtractor
//        EventExtractor extractor = EventExtractorFactory.create(config, ....);

        // TODO 3. Create ActionableEvents for all relevant entries.
    }
}
