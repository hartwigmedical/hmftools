package com.hartwig.hmftools.common.utils;


import java.util.List;

import com.google.common.collect.Lists;

public class GenericDataCollection {

    public static int GD_TYPE_DECIMAL = 0;
    public static int GD_TYPE_INTEGER = 1;
    public static int GD_TYPE_STRING = 2;

    private List<String> mFieldNames;

    private List<List<Integer>> mIntData;
    private List<List<Double>> mDecimalData;
    private List<List<String>> mStringData;
    private int mDataType;

    public GenericDataCollection(int dataType)
    {
        mFieldNames = Lists.newArrayList();

        mDataType = dataType;

        mDecimalData = null;
        mIntData = null;
        mStringData = null;

        if(dataType == GD_TYPE_DECIMAL)
            mDecimalData = Lists.newArrayList();
        else if(dataType == GD_TYPE_INTEGER)
            mIntData = Lists.newArrayList();
        else if(dataType == GD_TYPE_STRING)
            mStringData = Lists.newArrayList();
    }

    public int getDataType() { return mDataType; }

    public static boolean isValidType(int dataType)
    {
        return (dataType == GD_TYPE_DECIMAL || dataType == GD_TYPE_INTEGER || dataType == GD_TYPE_STRING);
    }

    public int getDataCount()
    {
        if(mDataType == GD_TYPE_DECIMAL)
            return mDecimalData.size();
        else if(mDataType == GD_TYPE_INTEGER)
            return mIntData.size();
        else if(mDataType == GD_TYPE_STRING)
            return mStringData.size();

        return 0;
    }

    public void setFieldNames(final List<String> fieldNames) { mFieldNames.addAll(fieldNames); }
    public final List<String> getFieldNames() { return mFieldNames; }

    public void addDecimalValues(final List<Double> dataValues) { mDecimalData.add(dataValues); }
    public void addIntValues(final List<Integer> dataValues) { mIntData.add(dataValues); }
    public void addStringValues(final List<String> dataValues) { mStringData.add(dataValues); }

    public final List<List<Double>> getData() { return mDecimalData; } // for backwards compatibility
    public final List<List<Double>> getDecimalData() { return mDecimalData; }
    public final List<List<Integer>> getIntData() { return mIntData; }
    public final List<List<String>> getStringData() { return mStringData; }

}
