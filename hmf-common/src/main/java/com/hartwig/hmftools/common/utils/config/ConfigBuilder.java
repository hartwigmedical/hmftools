package com.hartwig.hmftools.common.utils.config;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.FLAG;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.PATH;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.STRING;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ConfigBuilder
{
    private final String mConfigPrefix;
    private final List<ConfigItem> mItems;

    private static final String DEFAULT_CONFIG_PREFIX = "-";
    private static final Logger LOGGER = LogManager.getLogger(ConfigBuilder.class);

    public ConfigBuilder()
    {
        this(DEFAULT_CONFIG_PREFIX);
    }

    public ConfigBuilder(final String prefix)
    {
        mItems = Lists.newArrayList();
        mConfigPrefix = prefix;
    }

    public void addConfigItem(final ConfigItem item) { mItems.add(item); }

    public void addConfigItem(final String name, final String description, final String defaultValue)
    {
        addConfigItem(STRING, name, false, description, defaultValue);
    }

    public void addFlagItem(final String name, final String description)
    {
        addConfigItem(FLAG, name, false, description, null);
    }

    public void addConfigItem(final String name, final boolean required, final String description)
    {
        addConfigItem(STRING, name, required, description, null);
    }

    public void addPathItem(final String name, final boolean required, final String description)
    {
        addConfigItem(PATH, name, required, description, null);
    }

    public void addConfigItem(final String name, final boolean required, final String description, final String defaultValue)
    {
        addConfigItem(STRING, name, required, description, defaultValue);
    }

    public void addRequiredConfigItem(final ConfigItemType type, final String name, final String description)
    {
        addConfigItem(type, name, true, description, null);
    }

    public void addRequiredConfigItem(final String name, final String description)
    {
        addConfigItem(STRING, name, true, description, null);
    }

    public void addConfigItem(
            final ConfigItemType type, final String name, final boolean required, final String description, final String defaultValue)
    {
        mItems.add(new ConfigItem(type, name, required, description, defaultValue));
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

    public double getDecimal(final String name) { return getItem(name).decimal(); }
    public int getInteger(final String name) { return getItem(name).integer(); }
    public boolean hasFlag(final String name) { return getItem(name).bool(); }

    public boolean isValid()
    {
        boolean allValid = true;

        for(ConfigItem item : mItems)
        {
            if(item.missing())
            {
                LOGGER.error("missing config({}) desc({})", item.Name, item.Description);
                allValid = false;
            }
            else if(item.Type == PATH && !Files.exists(Paths.get(item.value())))
            {
                LOGGER.error("invalid config({}) path({})", item.Name, item.value());
                allValid = false;
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
                    LOGGER.error("invalid config({}:{}) value({})", item.Name, item.Type, item.value());
                    allValid = false;
                }
            }
        }

        return allValid;
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
                    LOGGER.error("config item({}) has invalid argument {}", matchedItem.Name, argument);
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
                    LOGGER.error("unexpected config item {}", argument);
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
                LOGGER.error("expected config value at invalid argument {}", argument);
                return false;
            }
        }

        return isValid();
    }

    public void logItems()
    {
        for(ConfigItem item : mItems)
        {
            LOGGER.info("{}: type({}) default({}) required({}) desc: {}",
                    item.Name, item.Type, item.defaultValue() != null ? item.defaultValue() : "none",
                    item.Required, item.Description);
        }
    }

    public void clearValues()
    {
        mItems.forEach(x -> x.clearValue());
    }
}
