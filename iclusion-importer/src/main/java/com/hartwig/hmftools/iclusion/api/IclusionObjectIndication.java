package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.squareup.moshi.Json;

class IclusionObjectIndication {

    @Json(name = "id")
    public String id;

    @Json(name = "indication_name")
    public String indicationName;

    @Json(name = "indication_name_full")
    public String indicationNameFull;

    @Json(name = "node_ids")
    public List<String> nodeIds;
}
