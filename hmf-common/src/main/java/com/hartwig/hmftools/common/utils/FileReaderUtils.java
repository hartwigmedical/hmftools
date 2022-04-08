package com.hartwig.hmftools.common.utils;

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
