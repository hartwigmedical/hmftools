package com.hartwig.hmftools.healthchecker.context;

import java.io.File;

import org.jetbrains.annotations.NotNull;

public final class CPCTRunContextFactory {

    private static final int CPCT_PATIENT_NAME_LENGTH = 12;

    private static final String CPCT_PATIENT_IDENTIFIER = "_CPCT";
    private static final String REF_SAMPLE_SUFFIX = "R";
    private static final String TUMOR_SAMPLE_SUFFIX = "T";

    private CPCTRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws MalformedRunDirException {
        final String fullRunName = removePath(runDirectory);

        final int patientPosition = fullRunName.indexOf(CPCT_PATIENT_IDENTIFIER) + 1;
        if ((patientPosition == 0) || (fullRunName.length() < (patientPosition + CPCT_PATIENT_NAME_LENGTH))) {
            throw new MalformedRunDirException(fullRunName);
        }

        final String patient = fullRunName.substring(patientPosition, patientPosition + CPCT_PATIENT_NAME_LENGTH);
        final String runName = fullRunName.substring(0, patientPosition + CPCT_PATIENT_NAME_LENGTH);

        final boolean hasPassedTests = fullRunName.length() == runName.length();

        final File[] runContents = new File(runDirectory).listFiles();
        assert runContents != null;

        String refSample = null;
        String tumorSample = null;
        for (final File content : runContents) {
            if (content.isDirectory() && content.getName().contains(patient)) {
                if (content.getName().contains(REF_SAMPLE_SUFFIX)) {
                    refSample = content.getName();
                } else if (content.getName().contains(TUMOR_SAMPLE_SUFFIX)) {
                    tumorSample = content.getName();
                }
            }
        }

        if (refSample == null || tumorSample == null) {
            throw new MalformedRunDirException(fullRunName);
        }

        return new RunContextImpl(runDirectory, runName, refSample, tumorSample, hasPassedTests);
    }

    @NotNull
    private static String removePath(@NotNull final String runDirectory) {
        String folderName = runDirectory;
        if (runDirectory.contains(File.separator)) {
            folderName = runDirectory.substring(runDirectory.lastIndexOf(File.separator) + 1, runDirectory.length());
        }
        return folderName;
    }
}
