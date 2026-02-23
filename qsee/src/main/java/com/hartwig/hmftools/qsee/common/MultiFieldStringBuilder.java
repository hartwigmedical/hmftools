package com.hartwig.hmftools.qsee.common;

import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import org.apache.commons.lang3.tuple.Pair;

public class MultiFieldStringBuilder
{
    private static final String FIELD_KEY_VALUE_SEPARATOR = "=";
    private static final String FIELD_SEPARATOR = ";";

    private final List<Pair<String, String>> mKeyValuePairs;

    public MultiFieldStringBuilder()
    {
        mKeyValuePairs = new ArrayList<>();
    }

    public MultiFieldStringBuilder add(String key, String value)
    {
        mKeyValuePairs.add(Pair.of(key, value));
        return this;
    }

    public String toString()
    {
        StringJoiner joiner = new StringJoiner(FIELD_SEPARATOR);

        for(Pair<String, String> pair : mKeyValuePairs)
        {
            String key = pair.getKey();
            String value = pair.getValue();

            if(value.isEmpty())
                continue;

            joiner.add(formSingleField(key, value));
        }

        return joiner.toString();
    }

    public static String formSingleField(String key, String value)
    {
        return key + FIELD_KEY_VALUE_SEPARATOR + value;
    }

    public static String formMultiField(String... keyValuePairs)
    {
        if(keyValuePairs.length % 2 != 0)
        {
            throw new IllegalArgumentException("Must provide an even number of arguments (key-value pairs)");
        }

        MultiFieldStringBuilder builder = new MultiFieldStringBuilder();
        for(int i = 0; i < keyValuePairs.length; i += 2)
        {
            String key = keyValuePairs[i];
            String value = keyValuePairs[i + 1];
            builder.add(key, value);
        }
        return builder.toString();
    }
}
