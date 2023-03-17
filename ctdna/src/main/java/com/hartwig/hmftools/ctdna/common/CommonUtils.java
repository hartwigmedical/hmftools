package com.hartwig.hmftools.ctdna.common;

import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class CommonUtils
{
    public static final Logger CT_LOGGER = LogManager.getLogger(CommonUtils.class);
    public static final String DELIMETER = ",";

    public static double medianIntegerValue(final List<Integer> values)
    {
        if(values.isEmpty())
            return 0;

        Collections.sort(values);

        int index = values.size() / 2;

        if((values.size() % 2) == 0)
        {
            return (values.get(index - 1) + values.get(index)) * 0.5;
        }
        else
        {
            return values.get(index);
        }
    }

}
