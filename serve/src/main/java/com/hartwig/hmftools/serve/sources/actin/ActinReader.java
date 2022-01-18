package com.hartwig.hmftools.serve.sources.actin;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.sources.actin.curation.ActinCurator;
import com.hartwig.hmftools.serve.sources.actin.filter.ActinFilter;
import com.hartwig.hmftools.serve.sources.actin.filter.ActinFilterEntry;
import com.hartwig.hmftools.serve.sources.actin.filter.ActinFilterFile;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActinReader {

    private static final Logger LOGGER = LogManager.getLogger(ActinReader.class);

    private ActinReader(){
    }

    @NotNull
    public static List<ActinEntry> read(@NotNull String actinTrialTsv, @NotNull String actinFilterTsv) throws IOException {
        LOGGER.info("Reading ACTIN trial database from {}", actinTrialTsv);
        List<ActinEntry> actinEntries = ActinFileReader.read(actinTrialTsv);
        LOGGER.info(" Read {} entries", actinEntries.size());

        LOGGER.info("Reading ACTIN filter entries from {}", actinFilterTsv);
        List<ActinFilterEntry> actinFilterEntries = ActinFilterFile.read(actinFilterTsv);
        LOGGER.info(" Read {} filter entries", actinFilterEntries.size());

        return filter(curate(actinEntries), actinFilterEntries);
    }

    @NotNull
    private static List<ActinEntry> curate(@NotNull List<ActinEntry> actinEntries) {
        ActinCurator curator = new ActinCurator();

        LOGGER.info("Curating {} ACTIN entries", actinEntries.size());
        List<ActinEntry> curatedEntries = curator.run(actinEntries);

        curator.reportUnusedCurationEntries();

        return curatedEntries;
    }

    @NotNull
    private static List<ActinEntry> filter(@NotNull List<ActinEntry> entries, @NotNull List<ActinFilterEntry> actinFilterEntries) {
        ActinFilter filter = new ActinFilter(actinFilterEntries);

        LOGGER.info("Filtering {} ACTIN entries", entries.size());
        List<ActinEntry> filteredEntries = filter.run(entries);
        LOGGER.info(" Finished ACTIN filtering. {} entries remaining, {} entries have been removed",
                filteredEntries.size(),
                entries.size() - filteredEntries.size());

        filter.reportUnusedFilterEntries();

        return filteredEntries;
    }
}
