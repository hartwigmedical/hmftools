package com.hartwig.hmftools.common.utils.config;

import java.util.ArrayList;
import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ConfigGroup
{
    private final String mName;
    private final String mDescription;
    private final List<ConfigItem> mItems;

    public ConfigGroup(@NotNull final String name, @NotNull final String description)
    {
        mItems = new ArrayList<>();
        mName = name;
        mDescription = description;
    }

    public void addConfigItem(@NotNull final ConfigItem item)
    {
        mItems.add(item);
    }

    public void addConfigItem(
            @NotNull final ConfigItemType type, @NotNull final String name, final boolean required, @NotNull final String description,
            @Nullable final String defaultValue)
    {
        addConfigItem(new ConfigItem(type, name, required, description, defaultValue));
    }

    @NotNull
    protected List<ConfigItem> getConfigItems()
    {
        return mItems;
    }
}
