package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;

import org.jetbrains.annotations.NotNull;

final class IclusionApiObjectMapper {

    private IclusionApiObjectMapper() {
    }

    @NotNull
    static List<IclusionTrial> fromApiObjects(@NotNull List<IclusionObjectStudy> studies, @NotNull List<IclusionObjectIndication> indications,
            @NotNull List<IclusionObjectGene> genes, @NotNull List<IclusionObjectVariant> variants) {
        return Lists.newArrayList();
    }
}
