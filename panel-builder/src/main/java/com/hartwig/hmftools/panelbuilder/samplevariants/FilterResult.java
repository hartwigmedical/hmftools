package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;

import org.jetbrains.annotations.Nullable;

// Generic monad for representing passing a filter or failing with a reason.
public record FilterResult(
        // Null if the filter is passed.
        @Nullable String failReason
)
{
    public static FilterResult pass()
    {
        return new FilterResult(null);
    }

    public static FilterResult fail(final String reason)
    {
        return new FilterResult(reason);
    }

    public static FilterResult condition(boolean cond, final String failReason)
    {
        return cond ? FilterResult.pass() : FilterResult.fail(failReason);
    }

    public boolean passed()
    {
        return failReason == null;
    }

    public <T> T map(final Supplier<T> onPassed, final Function<String, T> onFailed)
    {
        return passed() ? onPassed.get() : onFailed.apply(failReason);
    }

    public void unwrap(final Runnable onPassed, final Consumer<String> onFailed)
    {
        map(
                () ->
                {
                    onPassed.run();
                    return null;
                },
                reason ->
                {
                    onFailed.accept(reason);
                    return null;
                }
        );
    }

    public FilterResult thenFilter(final Supplier<FilterResult> filter)
    {
        return map(filter, failReason -> this);
    }

    public static FilterResult applyFilters(final List<Supplier<FilterResult>> filters)
    {
        return filters.stream().reduce(FilterResult::pass, (f1, f2) -> () -> f1.get().thenFilter(f2)).get();
    }
}
