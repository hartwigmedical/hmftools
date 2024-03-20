package com.hartwig.hmftools.peach;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

public class TestUtils
{
    @NotNull
    public static String getTestResourcePath(@NotNull String fileName)
    {
        return Resources.getResource(fileName).getPath();
    }
}
