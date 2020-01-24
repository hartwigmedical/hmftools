package com.hartwig.hmftools.iclusion.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;

import org.jetbrains.annotations.NotNull;

public final class IclusionTrialFile {

    private IclusionTrialFile() {
    }

    private static final String DELIMITER = "\t";

    public void write(@NotNull String iClusionTrialsTsv, @NotNull List<IclusionTrial> trials) throws IOException {
        Files.write(new File(iClusionTrialsTsv).toPath(), toLines(trials));
    }

    @NotNull
    private static List<String> toLines(@NotNull List<IclusionTrial> trials) {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        trials.stream().map(IclusionTrialFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("id")
                .add("acronym")
                .add("title")
                .add("eudra")
                .add("nct")
                .add("ipn")
                .add("ccmo")
                .add("tumorLocations")
                .add("mutations")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull IclusionTrial trial) {
        return new StringJoiner(DELIMITER).add(trial.id())
                .add(trial.acronym())
                .add(trial.title())
                .add(trial.eudra())
                .add(trial.nct())
                .add(trial.ipn())
                .add(trial.ccmo())
                .add("locations")
                .add("mutations")
                .toString();
    }
}
