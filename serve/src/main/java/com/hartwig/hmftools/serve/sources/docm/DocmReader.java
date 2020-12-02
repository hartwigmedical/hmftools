package com.hartwig.hmftools.serve.sources.docm;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.sources.docm.curation.DocmCurator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class DocmReader {

    private static final Logger LOGGER = LogManager.getLogger(DocmReader.class);

    private DocmReader() {
    }

    @NotNull
    public static List<DocmEntry> readAndCurate(@NotNull String docmTsv) throws IOException {
        LOGGER.info("Reading DoCM TSV from '{}'", docmTsv);
        List<DocmEntry> entries = DocmFileReader.read(docmTsv);
        LOGGER.info(" Read {} entries", entries.size());

        return curate(entries);
    }

    @NotNull
    private static List<DocmEntry> curate(@NotNull List<DocmEntry> entries) {
        DocmCurator curator = new DocmCurator();

        LOGGER.info("Curating {} DoCM entries", entries.size());
        List<DocmEntry> curatedEntries = curator.curate(entries);
        LOGGER.info(" Finished DoCM curation. {} entries remaining, {} entries have been removed",
                curatedEntries.size(),
                entries.size() - curatedEntries.size());

        curator.reportUnusedBlacklistEntries();

        return curatedEntries;
    }
}
