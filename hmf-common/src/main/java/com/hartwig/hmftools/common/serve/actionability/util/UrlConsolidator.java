package com.hartwig.hmftools.common.serve.actionability.util;

import java.util.Set;

import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public interface UrlConsolidator<T extends ActionableEvent> {

    @NotNull
    T stripUrls(@NotNull T instance);

    @NotNull
    T buildWithUrls(@NotNull T instance, @NotNull Set<String> urls);
}
