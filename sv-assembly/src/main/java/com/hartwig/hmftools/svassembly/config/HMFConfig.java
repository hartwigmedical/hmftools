package com.hartwig.hmftools.svassembly.config;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigItemType;
import com.hartwig.hmftools.svassembly.util.StringUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import java.io.File;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.file.Path;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;

public interface HMFConfig
{
    Logger LOGGER = LogManager.getLogger(HMFConfig.class);

    static void addConfig(final ConfigBuilder configBuilder, final Class<?> interfaceType)
    {
        for(final Method method : interfaceType.getDeclaredMethods())
        {
            if(Modifier.isStatic(method.getModifiers()))
                continue;

            @Nullable
            final CommandLine annotation = method.getAnnotation(CommandLine.class);
            final boolean isOptional = method.getAnnotation(Optional.class) != null;

            final Class<?> type = method.getReturnType();
            final String description = annotation == null ? "" : annotation.description();
            final String name =
                    annotation == null || annotation.name().isEmpty() ? StringUtils.toSnakeCase(method.getName()) : annotation.name();

            final ConfigItemType itemType;
            if(type == String.class || type.isEnum())
                itemType = ConfigItemType.STRING;
            else if(type == Integer.class || type == int.class)
                itemType = ConfigItemType.INTEGER;
            else if(type == Path.class || type == File.class)
                itemType = ConfigItemType.PATH;
            else if(type == Double.class || type == double.class || type == Float.class || type == float.class)
                itemType = ConfigItemType.DECIMAL;
            else if(type == Boolean.class || type == boolean.class)
                itemType = ConfigItemType.BOOL;
            else
                throw new IllegalArgumentException("Failed to config " + name + " with type " + type);

            @Nullable
            final String defaultValue = isOptional ? method.getAnnotation(Optional.class).defaultValue() : null;
            configBuilder.addConfigItem(itemType, name, !isOptional, description, defaultValue);
        }
    }

    @Nullable
    static <T, U> T load(final Map<String, Object> values, final Class<T> interfaceType, final U builder)
    {
        return load(values::containsKey, values::get, interfaceType, builder);
    }

    static <T, U> T load(final ConfigBuilder configBuilder, final Class<T> interfaceType, final U builder)
    {
        return load(configBuilder::hasValue, configBuilder::getValue, interfaceType, builder);
    }

    private static <T, U> T load(final Predicate<String> hasConfig, final Function<String, Object> configGetter,
            final Class<T> interfaceType, final U builder)
    {
        try
        {
            boolean hasErrors = false;
            for(final Method method : interfaceType.getDeclaredMethods())
            {
                if(Modifier.isStatic(method.getModifiers()))
                    continue;

                @Nullable
                final CommandLine annotation = method.getAnnotation(CommandLine.class);
                final boolean isOptional = method.getAnnotation(Optional.class) != null;
                final boolean mustBeExistingFile = method.getAnnotation(ExistingFile.class) != null;

                final Class<?> type = method.getReturnType();
                final Method setter = builder.getClass().getDeclaredMethod(method.getName(), type);

                final String name =
                        annotation == null || annotation.name().isEmpty() ? StringUtils.toSnakeCase(method.getName()) : annotation.name();

                if(!hasConfig.test(name))
                {
                    if(!isOptional)
                    {
                        hasErrors = true;
                        LOGGER.error("Configuration option {} is required", name);
                        continue;
                    }

                    final String rawDefault = method.getAnnotation(Optional.class).defaultValue();
                    final Object defaultValue = parseValue(name, rawDefault, type, mustBeExistingFile, true);
                    setter.invoke(builder, defaultValue);

                    continue;
                }

                final String rawValue = String.valueOf(configGetter.apply(name));
                @Nullable
                final Object parsedValue = parseValue(name, rawValue, type, mustBeExistingFile, false);
                if(parsedValue == null)
                    hasErrors = true;
                else
                    setter.invoke(builder, parsedValue);
            }

            //noinspection unchecked
            return hasErrors ? null : (T) builder.getClass().getDeclaredMethod("build").invoke(builder);
        }
        catch(final Exception exception)
        {
            throw new RuntimeException(exception);
        }
    }

    private static <T> Object parseValue(final String configName, final String value, final Class<T> type,
            final boolean mustBeExistingFile, final boolean throwErrors)
    {
        if (value.equals("null"))
            return null;
        else if(type == String.class)
            return value;
        else if(type.isEnum())
        {
            try
            {
                return type.getDeclaredMethod("from", String.class).invoke(null, value);
            }
            catch(final Exception firstAttempt)
            {
                try
                {
                    return type.getDeclaredMethod("valueOf", String.class).invoke(null, value);
                }
                catch(final Exception exception)
                {
                    throw new RuntimeException(exception);
                }
            }
        }
        else if(type == Integer.class || type == int.class)
        {
            try
            {
                return Integer.parseInt(value);
            }
            catch(final NumberFormatException exception)
            {
                if(throwErrors)
                    throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName, exception);
                LOGGER.error("Failed to parse value {} as integer for {}", value, configName);
            }
        }
        else if(type == File.class)
        {
            if(mustBeExistingFile)
            {
                if(!new File(value).exists())
                {
                    if(throwErrors)
                        throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName);
                    LOGGER.error("Could not find file {} supplied for {}", value, configName);
                    return null;
                }
            }

            return new File(value);
        }
        else if(type == Path.class)
        {
            if(mustBeExistingFile)
            {
                if(!new File(value).exists())
                {
                    if(throwErrors)
                        throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName + ", file not found.");
                    LOGGER.error("Could not find file {} supplied for {}", value, configName);
                    return null;
                }
            }

            return Path.of(value);
        }
        else if(type == Double.class || type == double.class)
        {
            try
            {
                return Double.parseDouble(value);
            }
            catch(final NumberFormatException exception)
            {
                if(throwErrors)
                    throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName);
                LOGGER.error("Failed to parse value {} as double for {}", value, configName);
                return null;
            }
        }
        else if(type == Float.class || type == float.class)
        {
            try
            {
                return Float.parseFloat(value);
            }
            catch(final NumberFormatException exception)
            {
                if(throwErrors)
                    throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName);
                LOGGER.error("Failed to parse value {} as double for {}", value, configName);
                return null;
            }
        }
        else if(type == Boolean.class || type == boolean.class)
        {
            try
            {
                return StringUtils.parseBoolean(value);
            }
            catch(final IllegalArgumentException exception)
            {
                if(throwErrors)
                    throw new IllegalArgumentException("Failed to parse value " + value + " for " + configName);
                LOGGER.error("Failed to parse value {} as boolean for {}", value, configName);
                return null;
            }
        }
        else
        {
            throw new IllegalArgumentException("Failed to config " + configName + " with type " + type);
        }
        return null;
    }
}
