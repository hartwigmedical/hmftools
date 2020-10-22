package com.hartwig.hmftools.healthchecker.runners;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.common.utils.io.path.PathPrefixSuffixFinder;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;
import com.hartwig.hmftools.healthchecker.result.ImmutableQCValue;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class FlagstatChecker implements HealthChecker {

    private static final Logger LOGGER = LogManager.getLogger(FlagstatFile.class);

    @NotNull
    private final String refSample;
    @NotNull
    private final String tumSample;
    @NotNull
    private final String flagstatDirectory;

    public FlagstatChecker(@NotNull final String refSample, @NotNull final String tumSample, @NotNull final String flagstatDirectory) {
        this.tumSample = tumSample;
        this.refSample = refSample;
        this.flagstatDirectory = flagstatDirectory;
    }

    public static String divideTwoStrings(String string1, String string2) {
        return String.valueOf(Double.valueOf(string1) / Double.valueOf(string2));
    }

    @NotNull
    @Override
    public List<QCValue> run() throws IOException {

        String refFlagstat = PathPrefixSuffixFinder.build().findPath(flagstatDirectory, refSample, ".flagstat").toString();
        String tumFlagstat = PathPrefixSuffixFinder.build().findPath(flagstatDirectory, tumSample, ".flagstat").toString();

        // Example flagstat lines
        // 323329219 + 0 in total (QC-passed reads + QC-failed reads)
        // 322517430 + 0 mapped (99.75%:N/A)

        List<String> ref_total_line = LineReader.build().readLines(new File(refFlagstat).toPath(), x -> x.contains("in total"));
        assert ref_total_line.size() == 1;
        String ref_total = ref_total_line.get(0).split(" ")[0];

        List<String> ref_mapped_line = LineReader.build().readLines(new File(refFlagstat).toPath(), x -> x.contains("mapped ("));
        assert ref_mapped_line.size() == 1;
        String ref_mapped = ref_mapped_line.get(0).split(" ")[0];

        List<String> tum_total_line = LineReader.build().readLines(new File(tumFlagstat).toPath(), x -> x.contains("in total"));
        assert tum_total_line.size() == 1;
        String tum_total = tum_total_line.get(0).split(" ")[0];

        List<String> tum_mapped_line = LineReader.build().readLines(new File(tumFlagstat).toPath(), x -> x.contains("mapped ("));
        assert tum_mapped_line.size() == 1;
        String tum_mapped = tum_mapped_line.get(0).split(" ")[0];

        String ref_prop = divideTwoStrings(ref_mapped, ref_total);
        String tum_prop = divideTwoStrings(tum_mapped, tum_total);

        return Lists.newArrayList(
            ImmutableQCValue.of(QCValueType.REF_PROPORTION_MAPPED, ref_prop),
            ImmutableQCValue.of(QCValueType.TUM_PROPORTION_MAPPED, tum_prop)
        );

    }
}
