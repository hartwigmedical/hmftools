package com.hartwig.hmftools.panelbuilder;

import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;

import org.jetbrains.annotations.Nullable;

// Generic monad for representing the result of a check/filter that: accepts; or rejects with a reason.
public record EvaluationResult<R>(
        // Null if accepted.
        @Nullable R rejectionInfo
)
{
    public static <R> EvaluationResult<R> accept()
    {
        return new EvaluationResult<>(null);
    }

    public static <R> EvaluationResult<R> reject(final R rejectionInfo)
    {
        return new EvaluationResult<>(rejectionInfo);
    }

    public static <R> EvaluationResult<R> condition(boolean cond, final R rejectionInfo)
    {
        return cond ? EvaluationResult.accept() : EvaluationResult.reject(rejectionInfo);
    }

    public boolean accepted()
    {
        return rejectionInfo == null;
    }

    public boolean rejected()
    {
        return rejectionInfo != null;
    }

    public <T> T map(final Supplier<T> onAccepted, final Function<R, T> onRejected)
    {
        return accepted() ? onAccepted.get() : onRejected.apply(rejectionInfo);
    }

    public void unwrap(final Runnable onAccepted, final Consumer<R> onRejected)
    {
        map(
                () ->
                {
                    onAccepted.run();
                    return null;
                },
                reason ->
                {
                    onRejected.accept(reason);
                    return null;
                }
        );
    }

    public EvaluationResult<R> then(final Supplier<EvaluationResult<R>> filter)
    {
        return map(filter, reason -> this);
    }

    public static <R> EvaluationResult<R> applyEvaluations(final List<Supplier<EvaluationResult<R>>> filters)
    {
        return filters.stream().reduce(EvaluationResult::accept, (f1, f2) -> () -> f1.get().then(f2)).get();
    }
}
