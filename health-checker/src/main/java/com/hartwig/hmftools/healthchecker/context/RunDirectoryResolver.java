package com.hartwig.hmftools.healthchecker.context;

import java.io.File;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class RunDirectoryResolver {

    private static final String RUN_NAME_SEPARATOR = "_";

    private static final String CPCT_PATIENT_IDENTIFIER = "_CPCT";
    private static final int CPCT_PATIENT_NAME_LENGTH = 12;
    private static final String CPCT_REF_SAMPLE_SUFFIX = "R";
    private static final String CPCT_TUMOR_SAMPLE_SUFFIX = "T";

    private static final String DRUP_PATIENT_IDENTIFIER = "_DRUP";
    private static final int DRUP_PATIENT_NAME_LENGTH = 12;
    private static final String DRUP_REF_SAMPLE_SUFFIX = "R";
    private static final String DRUP_TUMOR_SAMPLE_SUFFIX = "T";

    private static final String PMC_PATIENT_IDENTIFIER = "_PMC";
    private static final int PMC_PATIENT_NAME_LENGTH = 9;
    private static final String PMC_REF_SAMPLE_SUFFIX = "R";
    private static final String PMC_TUMOR_SAMPLE_SUFFIX = "T";

    private static final int SINGLE_SAMPLE_RUN_NAME_PARTS = 4;

    private RunDirectoryResolver() {
    }

    @NotNull
    static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        final String runName = removePath(runDirectory);

        if (isCPCTRun(runName)) {
            return somatic(runName, runDirectory, CPCT_PATIENT_IDENTIFIER, CPCT_PATIENT_NAME_LENGTH,
                    CPCT_REF_SAMPLE_SUFFIX, CPCT_TUMOR_SAMPLE_SUFFIX);
        } else if (isDRUPRun(runName)) {
            return somatic(runName, runDirectory, DRUP_PATIENT_IDENTIFIER, DRUP_PATIENT_NAME_LENGTH,
                    DRUP_REF_SAMPLE_SUFFIX, DRUP_TUMOR_SAMPLE_SUFFIX);
        } else if (isPMCRun(runName)) {
            return somatic(runName, runDirectory, PMC_PATIENT_IDENTIFIER, PMC_PATIENT_NAME_LENGTH,
                    PMC_REF_SAMPLE_SUFFIX, PMC_TUMOR_SAMPLE_SUFFIX);
        } else if (isSingleSample(runName)) {
            return singleSample(runName, runDirectory);
        } else {
            throw new HartwigException("Run name not supported by health checker: " + runName);
        }
    }

    private static boolean isCPCTRun(@NotNull final String runName) {
        return runName.contains(CPCT_PATIENT_IDENTIFIER);
    }

    private static boolean isDRUPRun(@NotNull final String runName) {
        return runName.contains(DRUP_PATIENT_IDENTIFIER);
    }

    private static boolean isPMCRun(@NotNull final String runName) {
        return runName.contains(PMC_PATIENT_IDENTIFIER);
    }

    private static boolean isSingleSample(@NotNull final String runName) {
        final String[] values = runName.split(RUN_NAME_SEPARATOR);
        if (values.length <= SINGLE_SAMPLE_RUN_NAME_PARTS) {
            return true;
        } else {
            final String identifier = RUN_NAME_SEPARATOR + values[values.length - 1];
            return !identifier.equals(CPCT_PATIENT_IDENTIFIER) && !identifier.equals(DRUP_PATIENT_IDENTIFIER);
        }
    }

    @NotNull
    private static RunContext singleSample(@NotNull final String runName, @NotNull final String runDirectory)
            throws MalformedRunDirException {
        final String[] parts = runName.split(RUN_NAME_SEPARATOR);
        final String finalRunName =
                parts[0] + RUN_NAME_SEPARATOR + parts[1] + RUN_NAME_SEPARATOR + parts[2] + RUN_NAME_SEPARATOR
                        + parts[3];

        final String sample = parts[SINGLE_SAMPLE_RUN_NAME_PARTS - 1];
        final File[] runContents = new File(runDirectory).listFiles();
        assert runContents != null;

        boolean sampleIsFound = false;
        for (final File content : runContents) {
            if (content.isDirectory() && content.getName().equals(sample)) {
                sampleIsFound = true;
            }
        }

        if (!sampleIsFound) {
            throw new MalformedRunDirException(runName);
        }

        final boolean hasPassedTests = parts.length == SINGLE_SAMPLE_RUN_NAME_PARTS;

        return new RunContextImpl(runDirectory, finalRunName, sample, Strings.EMPTY, hasPassedTests, false);
    }

    @NotNull
    private static RunContext somatic(@NotNull final String runName, @NotNull final String runDirectory,
            @NotNull final String patientIdentifier, final int patientNameLength,
            @NotNull final String refSampleSuffix, @NotNull final String tumorSampleSuffix)
            throws MalformedRunDirException {
        final int patientPosition = runName.indexOf(patientIdentifier) + 1;
        // KODU: somatic has to be called with the right identifier, so this should always be found!
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

        return new RunContextImpl(runDirectory, runNameUntilPatient, refSample, tumorSample, hasPassedTests, true);
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
