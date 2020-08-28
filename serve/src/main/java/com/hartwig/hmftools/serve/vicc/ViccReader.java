package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.serve.vicc.curation.CurationKey;
import com.hartwig.hmftools.serve.vicc.curation.FeatureCurator;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;
import com.hartwig.hmftools.vicc.selection.ViccQuerySelection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccReader.class);

    @NotNull
    public static List<ViccEntry> readAndCurate(@NotNull String viccJsonPath, @NotNull ViccQuerySelection querySelection)
            throws IOException {
        FeatureCurator curator = new FeatureCurator();

        LOGGER.info("Reading VICC json from '{}' with sources '{}'", viccJsonPath, querySelection.sourcesToFilterOn());
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJsonPath, querySelection);
        LOGGER.info(" Read {} entries. Starting curation", viccEntries.size());

        List<ViccEntry> curatedViccEntries = curator.curate(viccEntries);
        LOGGER.info(" Finished curation. {} curated entries remaining. {} entries have been removed.",
                curatedViccEntries.size(),
                viccEntries.size() - curatedViccEntries.size());

        LOGGER.info("Analyzing usage of curation configuration keys");
        Set<CurationKey> unusedCurationKeys = curator.unusedCurationKeys();
        for (CurationKey unusedKey : unusedCurationKeys) {
            LOGGER.warn(" Unused key found: '{}'", unusedKey);
        }

        LOGGER.info("Finished analyzing usage of curation configuration keys. Found {} unused configuration entries.",
                unusedCurationKeys.size());

        return curatedViccEntries;
    }
}
