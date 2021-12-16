package com.hartwig.hmftools.serve.sources.actin.reader;

import java.util.List;

import com.hartwig.hmftools.ckb.classification.CkbEventTypeExtractor;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.iclusion.classification.IclusionEventTypeExtractor;
import com.hartwig.hmftools.serve.sources.actin.ActinExtractor;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinClassificationConfig;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinEventTypeExtractor;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinEntry {

    @NotNull
    @Value.Derived
    public List<EventType> type() {
        return ActinEventTypeExtractor.classify(this);
    }

    @NotNull
    public abstract String trial();

    @NotNull
    public abstract ActinRule rule();

    @NotNull
    public abstract List<String> parameters();
}
