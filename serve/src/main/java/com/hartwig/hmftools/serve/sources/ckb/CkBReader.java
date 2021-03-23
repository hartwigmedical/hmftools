package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.sources.ckb.filter.CKBFilter;
import com.hartwig.hmftools.serve.sources.vicc.ViccReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CkBReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccReader.class);

    private CkBReader(){

    }

    @NotNull
    public static List<CkbEntry> filterRelevantEntries(@NotNull List<CkbEntry> ckbEntries) {

        return filter(ckbEntries);
    }

    @NotNull
    private static List<CkbEntry> filter(@NotNull List<CkbEntry> entries) {
        CKBFilter filter = new CKBFilter();

        LOGGER.info("Filtering {} CKB entries", entries.size());
        List<CkbEntry> filteredEntries = filter.run(entries);
        LOGGER.info(" Finished CKB filtering. {} entries remaining, {} entries have been removed",
                filteredEntries.size(),
                entries.size() - filteredEntries.size());

        filter.reportUnusedFilterEntries();

        return filteredEntries;
    }

}
