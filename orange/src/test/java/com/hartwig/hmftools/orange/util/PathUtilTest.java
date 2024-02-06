package com.hartwig.hmftools.orange.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.File;

import com.google.common.io.Resources;

import org.junit.Test;

public class PathUtilTest
{
    private static final String RUN_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String REAL_ALGO_DIRECTORY = RUN_DIRECTORY + File.separator + "purple";
    private static final String FAKE_ALGO_DIRECTORY = RUN_DIRECTORY + File.separator + "fake-algo";

    @Test
    public void shouldResolveOptionalPathsCorrectly()
    {
        assertNull(PathUtil.optionalPath(null));
        assertNull(PathUtil.optionalPath(FAKE_ALGO_DIRECTORY));
        assertEquals(REAL_ALGO_DIRECTORY, PathUtil.optionalPath(REAL_ALGO_DIRECTORY));
    }

    @Test
    public void shouldPassThroughActualPathsIfMandatory()
    {
        assertEquals(REAL_ALGO_DIRECTORY, PathUtil.mandatoryPath(REAL_ALGO_DIRECTORY));
    }

    @Test(expected = IllegalArgumentException.class)
    public void shouldThrowExceptionOnNonExistingMandatoryPath()
    {
        PathUtil.mandatoryPath(FAKE_ALGO_DIRECTORY);
    }
}