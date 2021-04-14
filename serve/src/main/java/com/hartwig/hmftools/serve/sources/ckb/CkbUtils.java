package com.hartwig.hmftools.serve.sources.ckb;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CkbUtils {

    private static final Logger LOGGER = LogManager.getLogger(CkbUtils.class);
    private static final String FIELD_DELIMITER = "\t";

    private CkbUtils() {
    }

    public static void printExtractionResults(@NotNull ExtractionResult result) {
        LOGGER.info("Analysis performed on CKB extraction result");
        LOGGER.info(" {} actionable hotspot records generated", result.actionableHotspots().size());
        LOGGER.info(" {} actionable range records generated", result.actionableRanges().size());
        LOGGER.info(" {} actionable gene records generated", result.actionableGenes().size());
        LOGGER.info(" {} actionable fusion records generated", result.actionableFusions().size());
        LOGGER.info(" {} actionable tumor characteristics records generated", result.actionableCharacteristics().size());
    }

    public static void writeEventsToTsv(@NotNull String eventTsv, @NotNull List<CkbEntry> entries) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("profileName").add("firstGene").add("firstVariant").add("type").toString();
        lines.add(header);

        for (CkbEntry entry : entries) {
            lines.add(new StringJoiner(FIELD_DELIMITER).add(entry.profileName())
                    .add(entry.variants().get(0).gene().geneSymbol())
                    .add(entry.variants().get(0).variant())
                    .add(entry.type().toString())
                    .toString());
        }

        LOGGER.info("Writing {} CKB features to {}", lines.size() - 1, eventTsv);
        Files.write(new File(eventTsv).toPath(), lines);
    }
}

