package com.hartwig.hmftools.healthchecker.context;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class RunDirectoryResolverTest {

    private static final String RESOURCE_DIR = "RunDirectoryResolver";

    private static final String VALID_CPCT_RUNDIR = "160101_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";
    private static final String VALID_CPCT_RUN_NAME = VALID_CPCT_RUNDIR;
    private static final String VALID_CPCT_REF = "CPCT12345678R";
    private static final String VALID_CPCT_TUMOR = "CPCT12345678TII";

    private static final String VALID_DRUP_RUNDIR = "170101_HMFregDRUP_FR10002000_FR10012001_DRUP01020005";
    private static final String VALID_DRUP_RUN_NAME = VALID_DRUP_RUNDIR;
    private static final String VALID_DRUP_REF = "DRUP01020005R";
    private static final String VALID_DRUP_TUMOR = "DRUP01020005T";

    private static final String VALID_PMC_RUNDIR = "170404_HMFregPMC_FR10002000_FR20003000_PMC010001";
    private static final String VALID_PMC_RUN_NAME = VALID_PMC_RUNDIR;
    private static final String VALID_PMC_REF = "PMC010001R";
    private static final String VALID_PMC_TUMOR = "PMC010001T";

    private static final String LOW_QUAL_CPCT_RUNDIR = "160102_HMFregCPCT_FR10002000_FR20003000_CPCT12345678_LowQual";
    private static final String LOW_QUAL_CPCT_RUN_NAME = "160102_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";
    private static final String LOW_QUAL_CPCT_REF = "CPCT12345678R";
    private static final String LOW_QUAL_CPCT_TUMOR = "CPCT12345678T";

    private static final String VALID_SINGLE_SAMPLE_RUNDIR = "170202_HMFreg0100_FR10002000_SAMPLE";
    private static final String VALID_SINGLE_SAMPLE_RUN_NAME = "170202_HMFreg0100_FR10002000_SAMPLE";
    private static final String VALID_SINGLE_SAMPLE_NAME = "SAMPLE";

    private static final String LOW_QUAL_SINGLE_SAMPLE_RUNDIR = "170303_HMFreg0100_FR10002000_SAMPLE_LowQual";
    private static final String LOW_QUAL_SINGLE_SAMPLE_RUN_NAME = "170303_HMFreg0100_FR10002000_SAMPLE";
    private static final String LOW_QUAL_SINGLE_SAMPLE_NAME = "SAMPLE";

    private static final String INVALID_CPCT_PATIENT_RUNDIR = "160103_HMFregCPCT_FR10002000_FR20003000_CPCT1234";
    private static final String MISSING_REF_CPCT_RUNDIR = "160104_HMFregCPCT_FR10002000_FR20003000_CPCT12345678";

    @Test
    public void resolveCorrectlyForValidCPCTRunWithTII() throws HartwigException {
        final RunContext runContext = RunDirectoryResolver.fromRunDirectory(toPath(VALID_CPCT_RUNDIR));
        assertEquals(VALID_CPCT_REF, runContext.refSample());
        assertEquals(VALID_CPCT_TUMOR, runContext.tumorSample());
        assertEquals(VALID_CPCT_RUN_NAME, runContext.setName());
        assertTrue(runContext.isSomaticRun());
    }

    @Test
    public void resolveCorrectlyForLowQualCPCTRun() throws HartwigException {
        final RunContext runContextLowQual = RunDirectoryResolver.fromRunDirectory(toPath(LOW_QUAL_CPCT_RUNDIR));
        assertEquals(LOW_QUAL_CPCT_REF, runContextLowQual.refSample());
        assertEquals(LOW_QUAL_CPCT_TUMOR, runContextLowQual.tumorSample());
        assertEquals(LOW_QUAL_CPCT_RUN_NAME, runContextLowQual.setName());
        assertTrue(runContextLowQual.isSomaticRun());
    }

    @Test
    public void worksForDRUP() throws HartwigException {
        final RunContext runContext = RunDirectoryResolver.fromRunDirectory(toPath(VALID_DRUP_RUNDIR));
        assertEquals(VALID_DRUP_REF, runContext.refSample());
        assertEquals(VALID_DRUP_TUMOR, runContext.tumorSample());
        assertEquals(VALID_DRUP_RUN_NAME, runContext.setName());
        assertTrue(runContext.isSomaticRun());
    }

    @Test
    public void worksForPMC() throws HartwigException {
        final RunContext runContext = RunDirectoryResolver.fromRunDirectory(toPath(VALID_PMC_RUNDIR));
        assertEquals(VALID_PMC_REF, runContext.refSample());
        assertEquals(VALID_PMC_TUMOR, runContext.tumorSample());
        assertEquals(VALID_PMC_RUN_NAME, runContext.setName());
        assertTrue(runContext.isSomaticRun());
    }

    @Test
    public void worksForValidSingleSample() throws HartwigException {
        final RunContext runContext = RunDirectoryResolver.fromRunDirectory(toPath(VALID_SINGLE_SAMPLE_RUNDIR));
        assertEquals(VALID_SINGLE_SAMPLE_NAME, runContext.refSample());
        assertEquals(VALID_SINGLE_SAMPLE_RUN_NAME, runContext.setName());
        assertFalse(runContext.isSomaticRun());
    }

    @Test
    public void worksForLowQualSingleSample() throws HartwigException {
        final RunContext runContext = RunDirectoryResolver.fromRunDirectory(toPath(LOW_QUAL_SINGLE_SAMPLE_RUNDIR));
        assertEquals(LOW_QUAL_SINGLE_SAMPLE_NAME, runContext.refSample());
        assertEquals(LOW_QUAL_SINGLE_SAMPLE_RUN_NAME, runContext.setName());
        assertFalse(runContext.isSomaticRun());
    }

    @Test(expected = MalformedRunDirException.class)
    public void exceptionOnCPCTRunDirWithTooShortPatientName() throws HartwigException {
        RunDirectoryResolver.fromRunDirectory(toPath(INVALID_CPCT_PATIENT_RUNDIR));
    }

    @Test(expected = MalformedRunDirException.class)
    public void exceptionOnCPCTRunDirWithMissingRef() throws HartwigException {
        RunDirectoryResolver.fromRunDirectory(toPath(MISSING_REF_CPCT_RUNDIR));
    }

    @NotNull
    private static String toPath(@NotNull final String runDirectory) {
        return Resources.getResource(RESOURCE_DIR + File.separator + runDirectory).getPath();
    }
}