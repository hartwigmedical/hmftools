package com.hartwig.hmftools.common.utils.file;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class FileReaderUtils
{
    public static Map<String,Integer> createFieldsIndexMap(final String fieldsHeader, final String delimiter)
    {
        final String[] items = fieldsHeader.split(delimiter,-1);
        final Map<String,Integer> fieldsIndexMap = Maps.newLinkedHashMap();

        for(int i = 0; i < items.length; ++i)
        {
            fieldsIndexMap.put(items[i], i);
        }

        return fieldsIndexMap;
    }

    public static String getValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final String defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return values[fieldsIndexMap.get(field)];
    }

    public static boolean getBoolValue(final Map<String,Integer> fieldsIndexMap, final String field, final String[] values)
    {
        return getBoolValue(fieldsIndexMap, field, false, values);
    }

    public static boolean getBoolValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final boolean defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return Boolean.parseBoolean(values[fieldsIndexMap.get(field)]);
    }

    public static int getIntValue(final Map<String,Integer> fieldsIndexMap, final String field, final String[] values)
    {
        return getIntValue(fieldsIndexMap, field, 0, values);
    }

    public static int getIntValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final int defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return Integer.parseInt(values[fieldsIndexMap.get(field)]);
    }

    public static long getLongValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final String[] values)
    {
        return getLongValue(fieldsIndexMap, field, (long)0, values);
    }

    public static long getLongValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final long defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return Long.parseLong(values[fieldsIndexMap.get(field)]);
    }

    public static double getDoubleValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final String[] values)
    {
        return getDoubleValue(fieldsIndexMap, field, 0.0, values);
    }

    public static double getDoubleValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final double defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return Double.parseDouble(values[fieldsIndexMap.get(field)]);
    }

    public static byte getByteValue(final Map<String,Integer> fieldsIndexMap, final String field, final String[] values)
    {
        return getByteValue(fieldsIndexMap, field, (byte)0, values);
    }

    public static byte getByteValue(
            final Map<String,Integer> fieldsIndexMap, final String field, final byte defaultValue, final String[] values)
    {
        if(!fieldsIndexMap.containsKey(field))
            return defaultValue;

        return Byte.parseByte(values[fieldsIndexMap.get(field)]);
    }

    public static String getValue(final List<String> lines, final String field, final String defaultValue, final String delimiter)
    {
        for(String line : lines)
        {
            String[] values = line.split(delimiter, -1);
            if(values[0].equals(field))
                return values[1];
        }

        return defaultValue;
    }

}
