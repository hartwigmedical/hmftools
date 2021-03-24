package com.hartwig.hmftools.serve.sources.ckb;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.serve.sources.ckb.filter.CKBFilter;
import com.hartwig.hmftools.serve.sources.iclusion.curation.IclusionCurator;
import com.hartwig.hmftools.serve.sources.iclusion.filter.IclusionFilter;
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

//    @NotNull
//    public static List<IclusionTrial> readAndCurate(@NotNull String iClusionTrialTsv) throws IOException {
//        LOGGER.info("Reading iClusion trial TSV from '{}'", iClusionTrialTsv);
//        List<IclusionTrial> trials = IclusionTrialFile.read(iClusionTrialTsv);
//        LOGGER.info(" Read {} trials", trials.size());
//
//        return filter(curate(trials));
//    }
//
//    @NotNull
//    private static List<CkbEntry> curate(@NotNull List<CkbEntry> trials) {
//        IclusionCurator curator = new IclusionCurator();
//
//        LOGGER.info("Curating {} iClusion trials", trials.size());
//        List<IclusionTrial> curatedTrials = curator.run(trials);
//        LOGGER.info(" Finished iClusion curation. {} trials remaining, {} trials have been removed",
//                curatedTrials.size(),
//                trials.size() - curatedTrials.size());
//
//        curator.reportUnusedCurationEntries();
//
//        return curatedTrials;
//    }

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
