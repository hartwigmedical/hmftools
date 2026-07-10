package com.hartwig.hmftools.compar.common.field;

public class UnsupportedFieldOverrideException extends RuntimeException
{
    public UnsupportedFieldOverrideException(final Field field, final String settingName)
    {
        super(String.format("field(%s) type(%s) does not support a %s override", field.name(), field.type(), settingName));
    }
}
