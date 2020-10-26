package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.squareup.moshi.Json;

public class IclusionObjectMutationCondition {
    @Json(name = "mutations")
    public List<IclusionObjectMutation> mutations;
    @Json(name = "logicType")
    public String logicType;
}
