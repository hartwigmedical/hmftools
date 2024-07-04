package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.PATHS;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.CONFIG_FILE_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.FLAG;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.PATH;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.STRING;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ConfigBuilder
{
    private final String mConfigPrefix;
    private final String mAppName;
    private final List<ConfigItem> mItems;
    private final Set<ErrorType> mErrors;
    private boolean mWarnOnRepeatedRegos;

    private static final String DEFAULT_CONFIG_PREFIX = "-";
    private static final String PRINT_HELP = "-help";
    private static final String PRINT_VERSION = "-version";

    private static final Logger LOGGER = LogManager.getLogger(ConfigBuilder.class);

    private enum ErrorType
    {
        INVALID_PATH,
        INVALID_TYPE,
        MISSING_REQUIRED,
        INCORRECT_ARGUMENT;
    }

    public ConfigBuilder()
    {
        this(DEFAULT_CONFIG_PREFIX, null);
    }

    public ConfigBuilder(final String appName)
    {
        this(DEFAULT_CONFIG_PREFIX, appName);
    }

    public ConfigBuilder(final String prefix, final String appName)
    {
        mAppName = appName;
        mItems = Lists.newArrayList();
        mErrors = Sets.newHashSet();
        mConfigPrefix = prefix;
        mWarnOnRepeatedRegos = true;
    }

    public void disableWarnOnRepeatedRegos() { mWarnOnRepeatedRegos = false; }

    public void addConfigItem(final ConfigItem item)
    {
        ConfigItem matched = getItem(item.Name, false);

        if(matched != null)
        {
            if(mWarnOnRepeatedRegos)
            {
                LOGGER.warn("registering config item({}) again", matched);
            }

            return;
        }

        mItems.add(item);
    }

    public void addConfigItem(
            final ConfigItemType type, final String name, final boolean required, final String description, final String defaultValue)
    {
        addConfigItem(new ConfigItem(type, name, required, description, defaultValue));
    }

    // convenience methods
    public void addConfigItem(final String name, final String description)
    {
        addConfigItem(STRING, name, false, description, null);
    }

    public void addRequiredDecimal(final String name, final String description)
    {
        addConfigItem(DECIMAL, name, true, description, null);
    }

    public void addDecimal(final String name, final String desc, double defaultValue)
    {
        addConfigItem(DECIMAL, name, false, format("%s, default=%.3g", desc, defaultValue), String.valueOf(defaultValue));
    }

    public void addRequiredInteger(final String name, final String description)
    {
        addConfigItem(INTEGER, name, true, description, null);
    }

    public void addInteger(final String name, final String description, int defaultValue)
    {
        addConfigItem(INTEGER, name, false, format("%s, default=%d", description, defaultValue), String.valueOf(defaultValue));
    }

    public void addFlag(final String name, final String description)
    {
        addConfigItem(FLAG, name, false, description, null);
    }

    public void addConfigItem(final String name, final boolean required, final String description)
    {
        addConfigItem(STRING, name, required, description, null);
    }

    public void addPath(final String name, final boolean required, final String description)
    {
        addConfigItem(PATH, name, required, description, null);
    }

    public void addPaths(final String name, final boolean required, final String description)
    {
        addConfigItem(PATHS, name, required, description, null);
    }

    // add a config path item with an optional path prefix set from another item
    public void addPrefixedPath(final String name, final boolean required, final String description, final String prefixConfig)
    {
        ConfigItem item = new ConfigItem(PATH, name, required, description, null);
        item.setPathPrefixName(prefixConfig);
        addConfigItem(item);
    }

    public void addConfigItem(final String name, final boolean required, final String description, final String defaultValue)
    {
        addConfigItem(STRING, name, required, description, defaultValue);
    }

    public void addRequiredConfigItem(final String name, final String description)
    {
        addConfigItem(STRING, name, true, description, null);
    }

    public ConfigItem getItem(final String name) { return getItem(name, true); }

    public ConfigItem getItem(final String name, boolean logWarning)
    {
        ConfigItem item = mItems.stream().filter(x -> x.Name.equals(name)).findFirst().orElse(null);

        if(item == null && logWarning)
        {
            LOGGER.warn("invalid config item({}) requested", name);
        }

        return item;
    }

    public String getValue(final String name) { return getItem(name).value(); }
    public boolean hasValue(final String name) { return getItem(name).hasValue(); }

    public boolean isRegistered(final String name) { return getItem(name, false) != null; }

    public String getValue(final String name, final String defaultValue)
    {
        // provides a default even if one hasn't been set
        ConfigItem item = getItem(name);
        return item.hasValue() ? item.value() : defaultValue;
    }

    public double getDecimal(final String name) { return getItem(name).decimal(); }
    public int getInteger(final String name) { return getItem(name).integer(); }
    public boolean hasFlag(final String name) { return getItem(name).bool(); }

    public static double getConfigDecimal(final ConfigBuilder configBuilder, final String configName, final double defaultValue)
    {
        if(configBuilder.hasValue(configName))
            return configBuilder.getDecimal(configName);

        return defaultValue;
    }

    public static int getConfigInteger(final ConfigBuilder configBuilder, final String configName, final int defaultValue)
    {
        if(configBuilder.hasValue(configName))
            return configBuilder.getInteger(configName);

        return defaultValue;
    }

    public boolean isValid()
    {
        for(ConfigItem item : mItems)
        {
            if(item.missing())
            {
                LOGGER.error("missing required config: {}: {}", item.Name, item.Description);
                mErrors.add(ErrorType.MISSING_REQUIRED);
            }
            else if(item.Type == PATH || item.Type == PATHS)
            {
                if(!isValidPath(item))
                {
                    LOGGER.error("invalid path for config: {} = {}", item.Name, item.value());
                    mErrors.add(ErrorType.INVALID_PATH);
                }
            }
            else
            {
                // check types
                try
                {
                    if(item.Type == DECIMAL)
                        item.decimal();
                    else if(item.Type == INTEGER)
                        item.integer();
                }
                catch(Exception e)
                {
                    LOGGER.error("invalid type for config: {}: {} = {}", item.Name, item.Type, item.value());
                    mErrors.add(ErrorType.INVALID_TYPE);
                }
            }
        }

        return mErrors.isEmpty();
    }

    private boolean isValidPath(final ConfigItem item)
    {
        if(!item.hasValue() || item.value().contains("*"))
            return true;

        List<String> paths = item.Type == PATHS ?
                Arrays.stream(item.value().split(CONFIG_FILE_DELIM, -1)).collect(Collectors.toList()) : Lists.newArrayList(item.value());

        ConfigItem pathPrefixConfigItem = item.pathPrefixName() != null ? getItem(item.pathPrefixName()) : null;

        for(String path : paths)
        {
            if(pathPrefixConfigItem != null && pathPrefixConfigItem.hasValue())
            {
                String pathPrefix = checkAddDirSeparator(pathPrefixConfigItem.value());
                path = pathPrefix + path;
            }

            if(!Files.exists(Paths.get(path)))
                return false;
        }

        return true;
    }

    public void checkAndParseCommandLine(final String[] args)
    {
        if(!parseCommandLine(args))
        {
            logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(this);

        printAppVersion(true);
    }

    public boolean parseCommandLine(final String[] args)
    {
        if(args == null)
            return false;

        if(args.length == 1 && checkHelpAndVersion(args[0]))
            return true;

        ConfigItem matchedItem = null;

        for(int i = 0; i < args.length; ++i)
        {
            String argument = args[i];

            if(matchedItem != null)
            {
                if(argument.startsWith(mConfigPrefix) && matchedItem.Type != DECIMAL && matchedItem.Type != INTEGER) // account for negative values
                {
                    LOGGER.error("config item({}) has invalid argument: {}", matchedItem.Name, argument);
                    mErrors.add(ErrorType.INCORRECT_ARGUMENT);
                    return false;
                }

                matchedItem.setValue(argument);
                matchedItem = null;
            }
            else if(argument.startsWith(mConfigPrefix))
            {
                String configName = argument.substring(mConfigPrefix.length());
                matchedItem = getItem(configName, false);

                if(matchedItem == null)
                {
                    LOGGER.error("unregistered config item: {}", argument);
                    mErrors.add(ErrorType.INCORRECT_ARGUMENT);
                    return false;
                }

                if(matchedItem.Type == FLAG)
                {
                    matchedItem.setValue(Boolean.TRUE.toString());
                    matchedItem = null;
                }
            }
            else
            {
                LOGGER.error("expecting config name but encountered invalid argument: {}", argument);
                mErrors.add(ErrorType.INCORRECT_ARGUMENT);
                return false;
            }
        }

        return isValid();
    }

    public void logInvalidDetails()
    {
        if(mErrors.isEmpty())
            return;

        if(mErrors.contains(ErrorType.INCORRECT_ARGUMENT))
        {
            logItems();
        }
    }

    public void logItems()
    {
        printItems(true);
    }

    private void printItems(boolean asLog)
    {
        List<String> output = Lists.newArrayList();

        for(ConfigItem item : mItems)
        {
            StringBuilder sb = new StringBuilder();
            sb.append(format(" -%s:", item.Name));

            if(item.Type != STRING)
                sb.append(format(" type(%s)", item.Type));

            if(item.defaultValue() != null && item.Type != FLAG)
                sb.append(format(" default(%s)", item.defaultValue()));

            if(item.Required)
                sb.append(" REQUIRED");

            sb.append(format(" desc(%s)", item.Description));

            output.add(sb.toString());
        }

        Collections.sort(output);

        output.add(0, "registered config items:");

        if(asLog)
            output.forEach(x -> LOGGER.info(x));
        else
            output.forEach(x -> System.out.println(x));
    }

    public void clearValues()
    {
        mErrors.clear();
        mItems.forEach(x -> x.clearValue());
    }

    private void printAppVersion(boolean asLog)
    {
        VersionInfo versionInfo = mAppName != null ? fromAppName(mAppName) : null;

        if(versionInfo != null)
        {
            String versionMessage = format("%s version %s", mAppName, versionInfo.version());
            if(asLog)
                LOGGER.info(versionMessage);
            else
                System.out.println(versionMessage);
        }
    }

    private boolean checkHelpAndVersion(final String argument)
    {
        if(!argument.equals(PRINT_HELP) && !argument.equals(PRINT_VERSION))
            return false;

        printAppVersion(false);

        if(argument.equals(PRINT_HELP))
        {
            printItems(false);
            System.exit(0);
        }
        else if(argument.equals(PRINT_VERSION))
        {
            if(mAppName == null)
            {
                LOGGER.warn("application name or version not set in config builder");
            }

            System.exit(0);
        }

        return false;
    }

    @VisibleForTesting
    public void setValue(final String name, final String value) { getItem(name, false).setValue(value); }
    public void setValue(final String name, final boolean value) { getItem(name, false).setValue(String.valueOf(value)); }
    public void setValue(final String name, final int value) { getItem(name, false).setValue(String.valueOf(value)); }
    public void setValue(final String name, final double value) { getItem(name, false).setValue(String.valueOf(value)); }
}
