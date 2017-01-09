package com.hartwig.hmftools.healthchecker.context;

import static org.junit.Assert.assertEquals;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ProductionRunContextFactoryTest {

    private static final String RESOURCE_DIR = "CPCTRunContextFactory";

    private static final String VALID_RUNDIR = "160101_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";
    private static final String VALID_RUN_NAME = VALID_RUNDIR;
    private static final String VALID_REF = "CPCT12345678R";
    private static final String VALID_TUMOR = "CPCT12345678TII";

    private static final String LOW_QUAL_RUNDIR = "160102_HMFregCPCT_FR10002000_FR20003000_CPCT12345678_LowQual";
    private static final String LOW_QUAL_RUN_NAME = "160102_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";
    private static final String LOW_QUAL_REF = "CPCT12345678R";
    private static final String LOW_QUAL_TUMOR = "CPCT12345678T";

    private static final String INVALID_PATIENT_RUNDIR = "160103_HMFregCPCT_FR10002000_FR20003000_CPCT1234";
    private static final String MISSING_REF_RUNDIR = "160104_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";

    @Test
    public void resolveCorrectlyForValidRunWithTII() throws HartwigException {
        final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(toPath(VALID_RUNDIR));
        assertEquals(VALID_REF, runContext.refSample());
        assertEquals(VALID_TUMOR, runContext.tumorSample());
        assertEquals(VALID_RUN_NAME, runContext.runName());
        assertEquals(true, runContext.hasPassedTests());
    }

    @Test
    public void resolveCorrectlyForLowQualRun() throws HartwigException {
        final RunContext runContextLowQual = ProductionRunContextFactory.fromRunDirectory(toPath(LOW_QUAL_RUNDIR));
        assertEquals(LOW_QUAL_REF, runContextLowQual.refSample());
        assertEquals(LOW_QUAL_TUMOR, runContextLowQual.tumorSample());
        assertEquals(LOW_QUAL_RUN_NAME, runContextLowQual.runName());
        assertEquals(false, runContextLowQual.hasPassedTests());
    }

    @Test(expected = MalformedRunDirException.class)
    public void exceptionOnRunDirWithTooShortPatientName() throws HartwigException {
        ProductionRunContextFactory.fromRunDirectory(toPath(INVALID_PATIENT_RUNDIR));
    }

    @Test(expected = MalformedRunDirException.class)
    public void exceptionOnRunDirWithMissingRef() throws HartwigException {
        ProductionRunContextFactory.fromRunDirectory(toPath(MISSING_REF_RUNDIR));
    }

    @NotNull
    private static String toPath(@NotNull final String runDirectory) {
        return Resources.getResource(RESOURCE_DIR + File.separator + runDirectory).getPath();
    }
}