package com.hartwig.hmftools.finding;

import java.nio.file.Path;

import com.hartwig.hmftools.finding.datamodel.VisualisationFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jspecify.annotations.Nullable;

public class VisualisationFileUtil
{
    private static final Logger LOGGER = LogManager.getLogger(VisualisationFileUtil.class);

    private VisualisationFileUtil()
    {
    }

    public static VisualisationFile create(final String path)
    {
        return new VisualisationFile(Path.of(path).getFileName().toString());
    }

    public static @Nullable VisualisationFile createNullable(@Nullable final String path)
    {
        if(path != null)
        {
            return create(path);
        }
        return null;
    }
}