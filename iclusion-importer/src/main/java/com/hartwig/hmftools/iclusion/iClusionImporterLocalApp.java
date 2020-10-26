package com.hartwig.hmftools.iclusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.iclusion.api.IclusionApiMain;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.qc.IclusionTrialChecker;

import org.jetbrains.annotations.NotNull;

public class iClusionImporterLocalApp {

    public static void main(@NotNull final String[] args) throws IOException {
        String endPoint = "https://iclusion.org/api";

        String config = System.getProperty("user.home") + "/hmf/tmp/iclusion";
        List<String> lines = Files.readAllLines(new File(config).toPath());

        IclusionCredentials credentials = ImmutableIclusionCredentials.builder()
                .clientId(lines.get(0))
                .clientSecret(lines.get(1))
                .username(lines.get(2))
                .password(lines.get(3))
                .build();

        List<IclusionTrial> trials = IclusionApiMain.readIclusionTrials(endPoint, credentials);

        IclusionTrialChecker.check(trials);
    }
}
