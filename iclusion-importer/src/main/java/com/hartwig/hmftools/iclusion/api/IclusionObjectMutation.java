package com.hartwig.hmftools.iclusion.api;

import com.squareup.moshi.Json;

class IclusionObjectMutation {

    @Json(name = "gene_id")
    public String geneId;

    @Json(name = "variant_id")
    public String variantId;

    @Json(name = "negation")
    public String negation;
}
