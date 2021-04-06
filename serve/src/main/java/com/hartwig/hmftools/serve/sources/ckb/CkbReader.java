package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.sources.ckb.curation.CkbCurator;
import com.hartwig.hmftools.serve.sources.ckb.filter.CkbFilter;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CkbReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccReader.class);

    private CkbReader(){

    }

    @NotNull
    public static List<CkbEntry> filterAndCurateRelevantEntries(@NotNull List<CkbEntry> ckbEntries, @Nullable Integer maxFilesToReadPerType) {

        return filter(curate(ckbEntries));
    }

    @NotNull
    public static List<CkbEntry> filterAndCurateRelevantEntries(@NotNull List<CkbEntry> ckbEntries) {

        return filter(curate(ckbEntries));
    }

    @NotNull
    private static List<CkbEntry> curate(@NotNull List<CkbEntry> ckbEntries) {
        CkbCurator curator = new CkbCurator();

        LOGGER.info("Curating {} CKB", ckbEntries.size());
        List<CkbEntry> curatedCKB = curator.run(ckbEntries);
        LOGGER.info(" Finished CKB curation. {} entries remaining, {} entries have been removed",
                ckbEntries.size(),
                ckbEntries.size() - curatedCKB.size());

        curator.reportUnusedCurationEntries();

        return curatedCKB;
    }

    @NotNull
    private static List<CkbEntry> filter(@NotNull List<CkbEntry> entries) {
        CkbFilter filter = new CkbFilter();

        LOGGER.info("Filtering {} CKB entries", entries.size());
        List<CkbEntry> filteredEntries = filter.run(entries);
        LOGGER.info(" Finished CKB filtering. {} entries remaining, {} entries have been removed",
                filteredEntries.size(),
                entries.size() - filteredEntries.size());

        filter.reportUnusedFilterEntries();

        return filteredEntries;
    }

}
