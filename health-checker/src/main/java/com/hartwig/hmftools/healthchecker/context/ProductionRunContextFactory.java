package com.hartwig.hmftools.healthchecker.context;

import java.io.File;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private static final String CPCT_PATIENT_IDENTIFIER = "_CPCT";
    private static final int CPCT_PATIENT_NAME_LENGTH = 12;
    private static final String CPCT_REF_SAMPLE_SUFFIX = "R";
    private static final String CPCT_TUMOR_SAMPLE_SUFFIX = "T";

    private static final String DRUP_PATIENT_IDENTIFIER = "_DRUP";
    private static final int DRUP_PATIENT_NAME_LENGTH = 12;
    private static final String DRUP_REF_SAMPLE_SUFFIX = "R";
    private static final String DRUP_TUMOR_SAMPLE_SUFFIX = "T";

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        final String runName = removePath(runDirectory);

        if (isCPCTRun(runName)) {
            return somatic(runName, runDirectory, CPCT_PATIENT_IDENTIFIER, CPCT_PATIENT_NAME_LENGTH,
                    CPCT_REF_SAMPLE_SUFFIX, CPCT_TUMOR_SAMPLE_SUFFIX);
        } else if (isDRUPRun(runName)) {
            return somatic(runName, runDirectory, DRUP_PATIENT_IDENTIFIER, DRUP_PATIENT_NAME_LENGTH,
                    DRUP_REF_SAMPLE_SUFFIX, DRUP_TUMOR_SAMPLE_SUFFIX);
        } else {
            throw new HartwigException("Currently only somatic pipelines from CPCT or DRUP are supported!");
        }
    }

    private static boolean isCPCTRun(@NotNull final String runName) {
        return runName.contains(CPCT_PATIENT_IDENTIFIER);
    }

    private static boolean isDRUPRun(@NotNull final String runName) {
        return runName.contains(DRUP_PATIENT_IDENTIFIER);
    }

    @NotNull
    private static RunContext somatic(@NotNull final String runName, @NotNull final String runDirectory,
            @NotNull final String patientIdentifier, final int patientNameLength,
            @NotNull final String refSampleSuffix, @NotNull final String tumorSampleSuffix)
            throws MalformedRunDirException {
        final int patientPosition = runName.indexOf(patientIdentifier) + 1;
        assert patientPosition > 0;

        if (runName.length() < (patientPosition + patientNameLength)) {
            throw new MalformedRunDirException(runName);
        }

        final String patient = runName.substring(patientPosition, patientPosition + patientNameLength);
        final String runNameUntilPatient = runName.substring(0, patientPosition + patientNameLength);

        final boolean hasPassedTests = runName.length() == runNameUntilPatient.length();

        final File[] runContents = new File(runDirectory).listFiles();
        assert runContents != null;

        String refSample = null;
        String tumorSample = null;
        for (final File content : runContents) {
            if (content.isDirectory() && content.getName().contains(patient)) {
                if (content.getName().substring(patientNameLength).contains(refSampleSuffix)) {
                    refSample = content.getName();
                } else if (content.getName().substring(patientNameLength).contains(tumorSampleSuffix)) {
                    tumorSample = content.getName();
                }
            }
        }

        if (refSample == null || tumorSample == null) {
            throw new MalformedRunDirException(runName);
        }

        return new RunContextImpl(runDirectory, runNameUntilPatient, refSample, tumorSample, hasPassedTests);
    }

    @NotNull
    private static String removePath(@NotNull final String runDirectory) {
        String folderName;
        if (runDirectory.contains(File.separator)) {
            folderName = runDirectory.substring(runDirectory.lastIndexOf(File.separator) + 1, runDirectory.length());
        } else {
            folderName = runDirectory;
        }
        return folderName;
    }
}
