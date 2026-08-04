package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;
import static com.hartwig.hmftools.compar.common.field.FieldInfo.findField;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.BoolFieldValue;
import com.hartwig.hmftools.compar.common.field.DoubleFieldValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.FieldValue;
import com.hartwig.hmftools.compar.common.field.IntFieldValue;
import com.hartwig.hmftools.compar.common.field.LongFieldValue;
import com.hartwig.hmftools.compar.common.field.StringFieldValue;

public abstract class ComparableItem
{
    protected final Map<String,FieldValue> mValues;

    public ComparableItem()
    {
        mValues = Maps.newHashMap();
    }

    public abstract CategoryType category();
    public abstract String key();

    public String geneName() { return ""; }
    public boolean reportable() { return true; }
    public boolean isPass() { return true; }

    public abstract boolean matches(final ComparableItem other);

    public Mismatch findMismatch(
            final ItemComparer comparer, final ComparableItem other, final MatchLevel matchLevel, final boolean includeMatches)
    {
        List<String> diffs = findDiffs(comparer, this, other);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    public Map<String,FieldValue> fieldValues() { return mValues; }

    public void addDoubleValue(final String fieldName, final Double value, final List<FieldInfo> fields)
    {
        mValues.put(fieldName, new DoubleFieldValue(findField(fieldName, fields), value));
    }

    public void addIntValue(final String fieldName, final Integer value, final List<FieldInfo> fields)
    {
        mValues.put(fieldName, new IntFieldValue(findField(fieldName, fields), value));
    }

    public void addLongValue(final String fieldName, final Long value, final List<FieldInfo> fields)
    {
        mValues.put(fieldName, new LongFieldValue(findField(fieldName, fields), value));
    }

    public void addBoolValue(final String fieldName, final boolean value, final List<FieldInfo> fields)
    {
        mValues.put(fieldName, new BoolFieldValue(findField(fieldName, fields), value));
    }

    public void addStringValue(final String fieldName, final String value, final List<FieldInfo> fields)
    {
        mValues.put(fieldName, new StringFieldValue(findField(fieldName, fields), value));
    }
}
