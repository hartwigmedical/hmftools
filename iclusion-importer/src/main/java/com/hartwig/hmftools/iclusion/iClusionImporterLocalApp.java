package com.hartwig.hmftools.iclusion;

import java.util.List;

import com.hartwig.hmftools.iclusion.api.IclusionApiMain;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.qc.IclusionTrialChecker;

import org.jetbrains.annotations.NotNull;

public class iClusionImporterLocalApp {

    public static void main(@NotNull final String[] args) {
        String endPoint = "https://iclusion.org/api";

        IclusionCredentials credentials = ImmutableIclusionCredentials.builder()
                .clientId("something")
                .clientSecret("something")
                .username("something")
                .password("something")
                .build();

        List<IclusionTrial> trials = IclusionApiMain.readIclusionTrials(endPoint, credentials);

        IclusionTrialChecker.check(trials);
    }
}
