package com.hartwig.hmftools.healthchecker.context;

import java.io.File;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;

public final class ProductionRunContextFactory {

    private static final int CPCT_PATIENT_NAME_LENGTH = 12;
    private static final String CPCT_PATIENT_IDENTIFIER = "_CPCT";
    private static final String CPCT_REF_SAMPLE_SUFFIX = "R";
    private static final String CPCT_TUMOR_SAMPLE_SUFFIX = "T";

    //    private static final int DRUP_PATIENT_NAME_LENGTH = 12;
    //    private static final String DRUP_PATIENT_IDENTIFIER = "_DRUP";
    //    private static final String DRUP_REF_SAMPLE_SUFFIX = "R";
    //    private static final String DRUP_TUMOR_SAMPLE_SUFFIX = "T";

    private ProductionRunContextFactory() {
    }

    @NotNull
    public static RunContext fromRunDirectory(@NotNull final String runDirectory) throws HartwigException {
        final String fullRunName = removePath(runDirectory);

        if (isCPCTRun(fullRunName)) {
            final int patientPosition = fullRunName.indexOf(CPCT_PATIENT_IDENTIFIER) + 1;
            if (fullRunName.length() < (patientPosition + CPCT_PATIENT_NAME_LENGTH)) {
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
                    if (content.getName().contains(CPCT_REF_SAMPLE_SUFFIX)) {
                        refSample = content.getName();
                    } else if (content.getName().contains(CPCT_TUMOR_SAMPLE_SUFFIX)) {
                        tumorSample = content.getName();
                    }
                }
            }

            if (refSample == null || tumorSample == null) {
                throw new MalformedRunDirException(fullRunName);
            }
            return new RunContextImpl(runDirectory, runName, refSample, tumorSample, hasPassedTests);
        } else {
            throw new HartwigException("CPCT is the only supported type of analysis for the health checker, sorry!");
        }
    }

    private static boolean isCPCTRun(@NotNull final String fullRunName) {
        return fullRunName.contains(CPCT_PATIENT_IDENTIFIER);
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
