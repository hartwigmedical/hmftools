package com.hartwig.hmftools.compar;

public class DiffData
{
    public final String FieldName;
    public final String RefValue;
    public final String OtherValue;

    public DiffData(final String fieldName, final String refValue, final String otherValue)
    {
        FieldName = fieldName;
        RefValue = refValue;
        OtherValue = otherValue;
    }
}
