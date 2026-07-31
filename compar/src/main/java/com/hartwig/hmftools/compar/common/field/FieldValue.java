package com.hartwig.hmftools.compar.common.field;

import static java.lang.String.format;

import java.util.List;

public abstract class FieldValue
{
    public final String Name;
    public final FieldCheck Check;
    public final String FormatString;

    public FieldValue(final String name, final FieldCheck check, final String formatString)
    {
        Name = name;
        Check = check;
        FormatString = formatString;
    }

    public abstract boolean checkDifference(final FieldValue other);

    public abstract String toString();
    public abstract String displayValue();

    public static void addDiffInfo(final FieldValue oldValue, final FieldValue newValue, final List<String> diffs)
    {
        diffs.add(format("%s(%s/%s)", oldValue.Name, oldValue.displayValue(), newValue.displayValue()));
    }
}
