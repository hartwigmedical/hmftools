package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.STRING;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigGroup;

import org.junit.Test;

public class ConfigBuilderTest
{
    @Test
    public void testConfigTypes()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(STRING, "item1", false, "desc1", "defaultValue");
        configBuilder.addConfigItem(DECIMAL, "item2", true, "desc2", "defaultValue");
        configBuilder.addConfigItem(INTEGER, "item3", false, "desc3", "defaultValue");
        configBuilder.addFlag("item4", "desc5");
        configBuilder.addRequiredConfigItem("item5", "desc5");

        String[] args = new String[9];
        int index = 0;
        args[index++] = "-item5";
        args[index++] = "val5";
        args[index++] = "-item4";
        args[index++] = "-item3";
        args[index++] = "10";
        args[index++] = "-item2";
        args[index++] = "10.56";
        args[index++] = "-item1";
        args[index++] = "val1";

        assertTrue(configBuilder.parseCommandLine(args));
        assertEquals("val1", configBuilder.getValue("item1"));
        assertEquals(10.56, configBuilder.getDecimal("item2"), 0.1);
        assertEquals(10, configBuilder.getInteger("item3"));
        assertTrue(configBuilder.hasFlag("item4"));
        assertEquals("val5", configBuilder.getValue("item5"));
    }

    @Test
    public void testInvalidConfig()
    {
        // missing config
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addRequiredConfigItem("item1", "desc1");
        configBuilder.addConfigItem("item2", false, "desc2", "defaultValue");

        String[] args = new String[2];
        int index = 0;
        args[index++] = "-item2";
        args[index++] = "val2";

        assertFalse(configBuilder.parseCommandLine(args));

        configBuilder.clearValues();

        // extra config
        args = new String[5];
        index = 0;
        args[index++] = "-item1";
        args[index++] = "val1";
        args[index++] = "-item2";
        args[index++] = "val2";
        args[index++] = "-test_flag";

        assertFalse(configBuilder.parseCommandLine(args));

        // invalid types
        configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(DECIMAL, "item1", true, "desc2", "defaultValue");
        configBuilder.addConfigItem(INTEGER, "item2", false, "desc2", "defaultValue");

        args = new String[4];
        index = 0;
        args[index++] = "-item1";
        args[index++] = "val";
        args[index++] = "-item2";
        args[index++] = "10.5";

        assertFalse(configBuilder.parseCommandLine(args));
    }

    @Test
    public void TestConfigGroup()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        // shared arg
        configBuilder.addConfigItem(STRING, "shared", true, "value for both groups", "item needed for both groups");

        // group of args handled together
        final ConfigGroup group1 = configBuilder.newConfigGroup("mode1", "group desc");
        group1.addConfigItem(STRING, "item1", true, "item desc", "default");

        // alternate group of args
        final ConfigGroup group2 = configBuilder.newConfigGroup("mode2", "group desc");
        group2.addConfigItem(STRING, "item2", true, "item desc", "default");

        // valid args for mode1. when groups are specified, a valid group must
        // be given as the first argument.
        String[] args1 = {
                "mode1",
                "-shared",
                "shared_value",
                "-item1",
                "val"
        };

        assertTrue(configBuilder.parseCommandLine(args1));

        // valid args for mode2
        String[] args2 = {
                "mode2",
                "-shared",
                "shared_value",
                "-item2",
                "val"
        };

        assertTrue(configBuilder.parseCommandLine(args2));
    }
}
