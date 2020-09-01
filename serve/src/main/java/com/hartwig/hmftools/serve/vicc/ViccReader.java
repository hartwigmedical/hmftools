package com.hartwig.hmftools.serve.vicc;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.serve.vicc.curation.FeatureCurator;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.reader.ViccJsonReader;
import com.hartwig.hmftools.vicc.selection.ImmutableViccQuerySelection;
import com.hartwig.hmftools.vicc.selection.ViccQuerySelection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViccReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccReader.class);

    @NotNull
    public static List<ViccEntry> readAndCurate(@NotNull String viccJson, @NotNull Set<ViccSource> sourcesToFilterOn, @Nullable Integer maxEntries)
            throws IOException {
        ViccQuerySelection querySelection =
                ImmutableViccQuerySelection.builder().sourcesToFilterOn(sourcesToFilterOn).maxEntriesToInclude(maxEntries).build();

        LOGGER.info("Reading VICC json from '{}' with sources '{}'", viccJson, querySelection.sourcesToFilterOn());
        List<ViccEntry> viccEntries = ViccJsonReader.readSelection(viccJson, querySelection);
        LOGGER.info(" Read {} entries", viccEntries.size());

        return curateViccEntries(viccEntries);
    }

    @NotNull
    public static List<ViccEntry> curateViccEntries(@NotNull List<ViccEntry> viccEntries) {
        FeatureCurator curator = new FeatureCurator();

        LOGGER.info("Curating {} entries.", viccEntries.size());
        List<ViccEntry> curatedViccEntries = curator.curate(viccEntries);
        LOGGER.info(" Finished curation. {} curated entries remaining. {} entries have been removed",
                curatedViccEntries.size(),
                viccEntries.size() - curatedViccEntries.size());

        curator.reportUnusedBlacklistEntries();

        return curatedViccEntries;
    }
}
