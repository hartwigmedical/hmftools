package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.List;

public abstract class FieldValue
{
    public final FieldInfo Field;

    public FieldValue(final FieldInfo field)
    {
        Field = field;
    }

    public String name() { return Field.Name; }
    public String formatString() { return Field.FormatString; }

    public abstract boolean hasDifference(final FieldValue other);

    public abstract String toString();
    public abstract String displayValue();

    public void addDiffInfo(final FieldValue oldValue, final FieldValue newValue, final List<String> diffs)
    {
        addDiffDisplayInfo(oldValue, newValue, diffs);
    }

    public static void addDiffDisplayInfo(final FieldValue oldValue, final FieldValue newValue, final List<String> diffs)
    {
        diffs.add(format("%s(%s/%s)", oldValue.name(), oldValue.displayValue(), newValue.displayValue()));
    }
}
