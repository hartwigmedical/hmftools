package com.hartwig.hmftools.svanalysis.annotators;

import java.util.List;

import com.google.common.collect.Lists;

public class ExternalSvData {

    private final String mId;
    private final List<String> mValues;

    public ExternalSvData(final String id, final List<String> values)
    {
        mId = id;
        mValues = values;
    }

    public final String getId() { return mId; }
    public final List<String> getValues() { return mValues; }

}
