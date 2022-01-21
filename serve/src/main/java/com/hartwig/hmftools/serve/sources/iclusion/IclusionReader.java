package com.hartwig.hmftools.serve.sources.iclusion;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.serve.sources.iclusion.curation.IclusionCurator;
import com.hartwig.hmftools.serve.sources.iclusion.filter.IclusionFilterEntry;
import com.hartwig.hmftools.serve.sources.iclusion.filter.IclusionFilterFile;
import com.hartwig.hmftools.serve.sources.iclusion.filter.IclusionFilter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class IclusionReader {

    private static final Logger LOGGER = LogManager.getLogger(IclusionReader.class);

    private IclusionReader() {
    }

    @NotNull
    public static List<IclusionTrial> readAndCurate(@NotNull String iClusionTrialTsv, @NotNull String iClusionFilterTsv) throws IOException {
        LOGGER.info("Reading iClusion trial TSV from '{}'", iClusionTrialTsv);
        List<IclusionTrial> trials = IclusionTrialFile.read(iClusionTrialTsv);
        LOGGER.info(" Read {} trials", trials.size());

        LOGGER.info("Reading iClusion filter entries from {}", iClusionFilterTsv);
        List<IclusionFilterEntry> iClusionFilterTrials = IclusionFilterFile.read(iClusionFilterTsv);
        LOGGER.info(" Read {} filter entries", iClusionFilterTrials.size());

        return filter(curate(trials), iClusionFilterTrials);
    }

    @NotNull
    private static List<IclusionTrial> curate(@NotNull List<IclusionTrial> trials) {
        IclusionCurator curator = new IclusionCurator();

        LOGGER.info("Curating {} iClusion trials", trials.size());
        List<IclusionTrial> curatedTrials = curator.run(trials);
        LOGGER.info(" Finished iClusion curation. {} trials remaining, {} trials have been removed",
                curatedTrials.size(),
                trials.size() - curatedTrials.size());

        curator.reportUnusedCurationEntries();

        return curatedTrials;
    }

    @NotNull
    private static List<IclusionTrial> filter(@NotNull List<IclusionTrial> trials, @NotNull List<IclusionFilterEntry> iClusionFilterTrials) {
        IclusionFilter filter = new IclusionFilter(iClusionFilterTrials);

        LOGGER.info("Filtering {} iClusion entries", trials.size());
        List<IclusionTrial> filteredTrials = filter.run(trials);
        LOGGER.info("  Finished iClusion filtering. {} trials remaining, {} entries have been removed",
                filteredTrials.size(),
                trials.size() - filteredTrials.size());

        filter.reportUnusedFilterEntries();

        return filteredTrials;
    }
}
