package com.hartwig.hmftools.svassembly;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.svassembly.config.HMFConfig;
import com.hartwig.hmftools.svassembly.processor.Processor;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.Nullable;

public class SVAssemblyApplication
{
    private static final Logger LOGGER = LogManager.getLogger(SVAssemblyApplication.class);

    public static void main(final String[] args)
    {
        final ConfigBuilder configBuilder = new ConfigBuilder();
        ConfigUtils.addLoggingOptions(configBuilder);
        SVAConfig.addConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        ConfigUtils.setLogLevel(configBuilder);
        Configurator.setLevel(LOGGER.getName(), Level.INFO);
        logVersion();

        @Nullable
        final SVAConfig config = HMFConfig.load(configBuilder, SVAConfig.class, ImmutableSVAConfig.builder());
        if(config == null)
            System.exit(1);

        LOGGER.info("Starting SVAssembly");

        final Context context = Context.create(config);

        new Processor(context).run();
    }

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("sv-assembly.version");
        LOGGER.info("SvAssembly version: {}", version.version());
    }
}
