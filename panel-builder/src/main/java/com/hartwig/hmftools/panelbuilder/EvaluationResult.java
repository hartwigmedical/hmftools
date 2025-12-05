package com.hartwig.hmftools.panelbuilder;

import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;

import org.jetbrains.annotations.Nullable;

// Generic monad for representing the result of a check/filter that: accepts; or rejects with a reason.
public record EvaluationResult(
        // Null if accepted.
        @Nullable String rejectionReason
)
{
    public EvaluationResult
    {
        if(rejectionReason != null && rejectionReason.isBlank())
        {
            throw new IllegalArgumentException("rejectionReason should not be blank");
        }
    }

    public static EvaluationResult accept()
    {
        return new EvaluationResult(null);
    }

    public static EvaluationResult reject(final String reason)
    {
        return new EvaluationResult(reason);
    }

    public static EvaluationResult condition(boolean cond, final String rejectionReason)
    {
        return cond ? EvaluationResult.accept() : EvaluationResult.reject(rejectionReason);
    }

    public boolean accepted()
    {
        return rejectionReason == null;
    }

    public boolean rejected()
    {
        return rejectionReason != null;
    }

    public <T> T map(final Supplier<T> onAccepted, final Function<String, T> onRejected)
    {
        return accepted() ? onAccepted.get() : onRejected.apply(rejectionReason);
    }

    public void unwrap(final Runnable onAccepted, final Consumer<String> onRejected)
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

    public EvaluationResult then(final Supplier<EvaluationResult> filter)
    {
        return map(filter, reason -> this);
    }

    public static EvaluationResult applyEvaluations(final List<Supplier<EvaluationResult>> filters)
    {
        return filters.stream().reduce(EvaluationResult::accept, (f1, f2) -> () -> f1.get().then(f2)).get();
    }
}
