package com.hartwig.hmftools.iclusion;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.iclusion.api.IclusionApiMain;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.iclusion.qc.IclusionTrialChecker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionImporter {

    private static final Logger LOGGER = LogManager.getLogger(IclusionImporter.class);

    @NotNull
    private final String apiEndPoint;
    @NotNull
    private final IclusionCredentials credentials;

    public IclusionImporter(@NotNull final String apiEndPoint, @NotNull final IclusionCredentials credentials) {
        this.apiEndPoint = apiEndPoint;
        this.credentials = credentials;
    }

    public void importToTsv(@NotNull String iClusionTrialTsv) throws IOException {
        List<IclusionTrial> trials = IclusionApiMain.readIclusionTrials(apiEndPoint, credentials);

        IclusionTrialChecker.check(trials);

        IclusionTrialFile.write(iClusionTrialTsv, trials);

        LOGGER.info("iClusion Importer has written {} trials to {}", trials.size(), iClusionTrialTsv);
    }
}
