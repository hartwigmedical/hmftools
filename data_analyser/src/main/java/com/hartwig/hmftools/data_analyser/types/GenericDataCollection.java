package com.hartwig.hmftools.data_analyser.types;

import java.util.List;

import com.google.common.collect.Lists;

public class GenericDataCollection {

    private List<String> mFieldNames;

    private List<List<Double>> mData;

    public GenericDataCollection()
    {
        mFieldNames = Lists.newArrayList();
        mData = Lists.newArrayList();
    }

    public void setFieldNames(final List<String> fieldNames) { mFieldNames.addAll(fieldNames); }
    public final List<String> getFieldNames() { return mFieldNames; }

    public void addDataValues(final List<Double> dataValues) { mData.add(dataValues); }
    public final List<List<Double>> getData() { return mData; }

}
