package com.hartwig.hmftools.common.sv;

import static java.lang.String.format;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantImpl implements StructuralVariant
{
    public String toString()
    {
        String coords = format("%s:%d:%d", start().chromosome(), start().position(), start().orientation());

        if(end() != null)
            coords += format(" - %s:%d:%d", end().chromosome(), end().position(), end().orientation());

        return format("id(%s) type(%s) coords(%s) filter(%s) qual(%.0f)",
                id(), type(), coords, filter(), qualityScore());
    }

}
