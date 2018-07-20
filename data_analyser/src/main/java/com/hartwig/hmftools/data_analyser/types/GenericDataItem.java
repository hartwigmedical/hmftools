package com.hartwig.hmftools.data_analyser.types;

import java.util.List;

import com.google.common.collect.Lists;

public class GenericDataItem {

    private int mFieldCount;
    private List<String> mFieldNames;

    public GenericDataItem()
    {
        mFieldCount = 0;
        mFieldNames = Lists.newArrayList();
    }

}
