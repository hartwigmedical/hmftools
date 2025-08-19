package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.List;
import java.util.function.Supplier;

import org.jetbrains.annotations.Nullable;

// Generic monad for representing passing a filter or failing with a reason.
public record FilterResult(
        @Nullable String failReason
)
{
    public static FilterResult pass()
    {
        return new FilterResult(null);
    }

    public static FilterResult fail(String reason)
    {
        return new FilterResult(reason);
    }

    public static FilterResult condition(boolean cond, String failReason)
    {
        if(cond)
        {
            return FilterResult.pass();
        }
        else
        {
            return FilterResult.fail(failReason);
        }
    }

    public boolean passed()
    {
        return failReason == null;
    }

    public FilterResult thenFilter(final Supplier<FilterResult> filter)
    {
        return passed() ? filter.get() : this;
    }

    public static FilterResult applyFilters(final List<Supplier<FilterResult>> filters)
    {
        FilterResult result = FilterResult.pass();
        for(Supplier<FilterResult> filter : filters)
        {
            result = result.thenFilter(filter);
        }
        return result;
    }
}
