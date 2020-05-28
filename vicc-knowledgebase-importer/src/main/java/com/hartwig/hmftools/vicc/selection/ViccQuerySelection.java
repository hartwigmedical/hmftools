package com.hartwig.hmftools.vicc.selection;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccQuerySelection {

    // If null: Include all sources
    @Nullable
    public abstract List<ViccSource> sourcesToFilterOn();

    // If null: Include all entries
    @Nullable
    public abstract Integer maxEntriesToInclude();
}
