package com.hartwig.hmftools.linx;

import java.io.BufferedWriter;

public interface CohortFileInterface
{
    String fileType();
    BufferedWriter createWriter(final String outputDir);
}
