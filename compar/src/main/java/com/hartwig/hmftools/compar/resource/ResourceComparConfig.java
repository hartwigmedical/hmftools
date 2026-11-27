package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ResourceComparConfig
{
    public final String OldResourcesDir;
    public final String NewResourcesDir;
    public final String OutputDir;
    public final String OutputId;

    public static final String OLD_RESOURCES_DIR = "old_resources_dir";
    public static final String NEW_RESOURCES_DIR = "new_resources_dir";

    public ResourceComparConfig(final ConfigBuilder configBuilder)
    {
        OldResourcesDir = configBuilder.getValue(OLD_RESOURCES_DIR);
        NewResourcesDir = configBuilder.getValue(NEW_RESOURCES_DIR);
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(OLD_RESOURCES_DIR, true, "Old resources directory, searched with its subdirectories");
        configBuilder.addPath(NEW_RESOURCES_DIR, true, "New resources directory, searched with its subdirectories");
        addOutputOptions(configBuilder, true);
        addLoggingOptions(configBuilder);
    }
}
