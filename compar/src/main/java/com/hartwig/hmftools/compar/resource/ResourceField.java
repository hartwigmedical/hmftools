package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;

import com.hartwig.hmftools.compar.common.field.BoolFieldValue;
import com.hartwig.hmftools.compar.common.field.DoubleFieldValue;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.FieldValue;
import com.hartwig.hmftools.compar.common.field.StringFieldValue;

import org.jetbrains.annotations.Nullable;

public class ResourceField
{
    public final FieldInfo Info;
    public final ValueType Type;
    @Nullable
    public final String DefaultValue;

    public enum ValueType
    {
        STRING,
        BOOLEAN,
        DOUBLE,
        GERMLINE_REPORTING
    }

    private static final String DISPLAY_FORMAT = "%s";

    private ResourceField(final FieldInfo info, final ValueType type, @Nullable final String defaultValue)
    {
        Info = info;
        Type = type;
        DefaultValue = defaultValue;
    }

    public static ResourceField compared(final String name, final ValueType type)
    {
        return new ResourceField(new FieldInfo(name, new FieldCheck(true), DISPLAY_FORMAT), type, null);
    }

    public static ResourceField compared(final String name, final ValueType type, final String defaultValue)
    {
        return new ResourceField(new FieldInfo(name, new FieldCheck(true), DISPLAY_FORMAT), type, defaultValue);
    }

    public static ResourceField displayOnly(final String name)
    {
        return new ResourceField(FieldInfo.displayOnly(name, DISPLAY_FORMAT), ValueType.STRING, null);
    }

    public String name() { return Info.Name; }

    public FieldValue createValue(final String value)
    {
        return switch(Type)
        {
            case STRING -> new StringFieldValue(Info, value);
            case BOOLEAN -> new BoolFieldValue(Info, Boolean.parseBoolean(value));
            case DOUBLE -> new DoubleFieldValue(Info, Double.parseDouble(value));
            case GERMLINE_REPORTING -> new StringFieldValue(Info, germlineReporting(value));
        };
    }

    private static String germlineReporting(final String value)
    {
        if(value.equalsIgnoreCase(Boolean.TRUE.toString()))
            return ANY.toString();

        if(value.equalsIgnoreCase(Boolean.FALSE.toString()))
            return NONE.toString();

        return value.toUpperCase();
    }
}
