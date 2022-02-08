package com.hartwig.hmftools.common.utils.config;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.Level;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;
import com.hartwig.hmftools.common.utils.ConfigUtils;

//
// This is for use with jcommander:
//
// @ParametersDelegate
// private final LoggingOptions mLoggingOptions = new LoggingOptions();
//
public class LoggingOptions
{
    // we need to define a converter for log4j level
    static class Log4jLevelTypeConverter implements IStringConverter<Level>
    {
        @Override
        public Level convert(String value)
        {
            return Level.valueOf(value);
        }
    }

    @Parameter(names = "-" + ConfigUtils.LOG_DEBUG, description = "Log verbose")
    public boolean LogDebug;

    @Parameter(names = "-" + ConfigUtils.LOG_LEVEL,
               description = "Specify log level (values: [OFF, FATAL, ERROR, WARN, INFO, DEBUG, TRACE, ALL])",
               converter = Log4jLevelTypeConverter.class)
    public org.apache.logging.log4j.Level LogLevel = null;

    public void setLogLevel()
    {
        if (LogDebug)
        {
            Configurator.setRootLevel(org.apache.logging.log4j.Level.DEBUG);
        }
        else if(LogLevel != null)
        {
            Configurator.setRootLevel(LogLevel);
        }
    }
}
