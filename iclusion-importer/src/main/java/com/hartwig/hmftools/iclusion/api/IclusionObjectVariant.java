package com.hartwig.hmftools.iclusion.api;

import com.squareup.moshi.Json;

class IclusionObjectVariant {

    @Json(name = "id")
    public String id;

    @Json(name = "variant_name")
    public String variantName;
}
