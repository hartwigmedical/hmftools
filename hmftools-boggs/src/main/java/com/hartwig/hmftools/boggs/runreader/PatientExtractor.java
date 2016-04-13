package com.hartwig.hmftools.boggs.runreader;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.boggs.PatientData;
import com.hartwig.hmftools.boggs.SampleData;
import com.hartwig.hmftools.boggs.flagstatreader.Flagstat;
import com.hartwig.hmftools.boggs.flagstatreader.FlagstatParser;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Collection;

public class PatientExtractor {

    private static final String SAMPLE_PREFIX = "CPCT";
    private static final String REF_SAMPLE_SUFFIX = "R";
    private static final String TUMOR_SAMPLE_SUFFIX = "T";
    private static final String FLAGSTAT_SUFFIX = ".flagstat";

    @NotNull
    private final FlagstatParser flagstatParser;

    public PatientExtractor(@NotNull FlagstatParser flagstatParser) {
        this.flagstatParser = flagstatParser;
    }

    @NotNull
    public PatientData extractFromRunDirectory(@NotNull String runDirectory) throws IOException {
        File directory = new File(runDirectory);

        SampleData refSample = extractSample(directory, refSampleFilter());
        SampleData tumorSample = extractSample(directory, tumorSampleFilter());

        return new PatientData(refSample, tumorSample);
    }

    @NotNull
    private SampleData extractSample(@NotNull File runDirectory, @NotNull FilenameFilter filter)
            throws IOException {
        File[] samples = runDirectory.listFiles(filter);

        assert samples.length == 1;

        String externalID = samples[0].getName();

        File flagstatDir = new File(samples[0].getPath() + File.separator + "mapping" + File.separator);
        File[] flagstats = flagstatDir.listFiles(flagstatFilter());

        Collection<Flagstat> mappingFlagstats = Lists.newArrayList();
        Flagstat markdupFlagstat = null;
        Flagstat realignedFlagstat = null;
        for (File flagstat : flagstats) {
            Flagstat parsedFlagstat = flagstatParser.parse(flagstat);
            String name = flagstat.getName();
            if (name.contains("dedup.realigned")) {
                realignedFlagstat = parsedFlagstat;
            } else if (name.contains("dedup")) {
                markdupFlagstat = parsedFlagstat;
            } else {
                mappingFlagstats.add(parsedFlagstat);
            }
        }

        assert markdupFlagstat != null;
        assert realignedFlagstat != null;

        return new SampleData(externalID, mappingFlagstats, markdupFlagstat, realignedFlagstat);
    }

    @NotNull
    private static FilenameFilter refSampleFilter() {
        return new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.startsWith(SAMPLE_PREFIX) && name.endsWith(REF_SAMPLE_SUFFIX);
            }
        };
    }

    @NotNull
    private static FilenameFilter tumorSampleFilter() {
        return new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.startsWith(SAMPLE_PREFIX) && name.endsWith(TUMOR_SAMPLE_SUFFIX);
            }
        };
    }

    @NotNull
    private static FilenameFilter flagstatFilter() {
        return new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(FLAGSTAT_SUFFIX);
            }
        };
    }
}
