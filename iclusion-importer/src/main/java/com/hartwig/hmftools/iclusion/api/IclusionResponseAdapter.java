package com.hartwig.hmftools.iclusion.api;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.squareup.moshi.FromJson;

import org.jetbrains.annotations.NotNull;

public class IclusionResponseAdapter {

    @FromJson
    @NotNull
    List<IclusionIndication> indications(@NotNull Map<String, IclusionIndication> json) {
        System.out.println("converting indications");
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionGene> genes(@NotNull Map<String, IclusionGene> json) {
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionVariant> variants(@NotNull Map<String, IclusionVariant> json) {
        return Lists.newArrayList(json.values());
    }

    @FromJson
    @NotNull
    List<IclusionStudy> studies(@NotNull Map<String, IclusionStudy> json) {
        return Lists.newArrayList(json.values());
    }
}
