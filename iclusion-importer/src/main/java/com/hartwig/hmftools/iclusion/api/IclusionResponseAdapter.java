package com.hartwig.hmftools.iclusion.api;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.squareup.moshi.FromJson;

import org.jetbrains.annotations.NotNull;

class IclusionResponseAdapter {

    @FromJson
    @NotNull
    List<IclusionObjectIndication> indications(@NotNull Map<String, IclusionObjectIndication> json) {
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionObjectGene> genes(@NotNull Map<String, IclusionObjectGene> json) {
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionObjectVariant> variants(@NotNull Map<String, IclusionObjectVariant> json) {
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionObjectStudy> studies(@NotNull Map<String, IclusionObjectStudy> json) {
        return Lists.newArrayList(json.values());
    }
}
