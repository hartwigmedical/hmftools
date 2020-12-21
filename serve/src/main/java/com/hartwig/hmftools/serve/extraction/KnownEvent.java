package com.hartwig.hmftools.serve.extraction;

import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;

public interface KnownEvent {

    @NotNull
    Set<Knowledgebase> sources();

}
