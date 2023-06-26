package com.hartwig.hmftools.common.utils.config;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.FLAG;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.PATH;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.STRING;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ConfigBuilder
{
    private final String mConfigPrefix;
    private final List<ConfigItem> mItems;
    private final Set<ErrorType> mErrors;

    private static final String DEFAULT_CONFIG_PREFIX = "-";
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
        this(DEFAULT_CONFIG_PREFIX);
    }

    public ConfigBuilder(final String prefix)
    {
        mItems = Lists.newArrayList();
        mErrors = Sets.newHashSet();
        mConfigPrefix = prefix;
    }

    public void addConfigItem(final ConfigItem item)
    {
        ConfigItem matched = mItems.stream().filter(x -> x.Name.equals(item.Name)).findFirst().orElse(null);

        if(matched != null)
        {
            LOGGER.warn("registering config item({}) again", matched);
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

    public void addConfigItem(final String name, final boolean required, final String description, final String defaultValue)
    {
        addConfigItem(STRING, name, required, description, defaultValue);
    }

    public void addRequiredConfigItem(final String name, final String description)
    {
        addConfigItem(STRING, name, true, description, null);
    }

    public ConfigItem getItem(final String name)
    {
        ConfigItem item = mItems.stream().filter(x -> x.Name.equals(name)).findFirst().orElse(null);

        if(item == null)
        {
            LOGGER.warn("invalid config item({}) requested", name);
        }

        return item;
    }

    public String getValue(final String name) { return getItem(name).value(); }
    public boolean hasValue(final String name) { return getItem(name).hasValue(); }

    public String getValue(final String name, final String defaultValue)
    {
        // provides a default even if one hasn't been set
        ConfigItem item = getItem(name);
        return item.hasValue() ? item.value() : defaultValue;
    }

    public double getDecimal(final String name) { return getItem(name).decimal(); }
    public int getInteger(final String name) { return getItem(name).integer(); }
    public boolean hasFlag(final String name) { return getItem(name).bool(); }

    public boolean isValid()
    {
        for(ConfigItem item : mItems)
        {
            if(item.missing())
            {
                LOGGER.error("missing required config: {}: {}", item.Name, item.Description);
                mErrors.add(ErrorType.MISSING_REQUIRED);
            }
            else if(item.Type == PATH)
            {
                if(item.hasValue() && !item.value().contains("*") && !Files.exists(Paths.get(item.value())))
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

    public boolean parseCommandLine(final String[] args)
    {
        if(args == null)
            return false;

        ConfigItem matchedItem = null;

        for(int i = 0; i < args.length; ++i)
        {
            String argument = args[i];

            if(matchedItem != null)
            {
                if(argument.startsWith(mConfigPrefix))
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
                matchedItem = mItems.stream().filter(x -> x.Name.equals(configName)).findFirst().orElse(null);

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
        LOGGER.info("{} registered config items:", mItems.size());

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

            LOGGER.info(sb.toString());
        }
    }

    public void clearValues()
    {
        mErrors.clear();
        mItems.forEach(x -> x.clearValue());
    }
}
