package com.hartwig.hmftools.esvee.util;

import java.util.LinkedHashMap;
import java.util.Map;

public class MapUtils
{
    public static Map<String, Object> mapOf(final Object... params)
    {
        final Map<String, Object> result = new LinkedHashMap<>();
        for(int i = 0; i < params.length; i += 2)
            result.put(params[i].toString(), params[i + 1]);
        return result;
    }
}
