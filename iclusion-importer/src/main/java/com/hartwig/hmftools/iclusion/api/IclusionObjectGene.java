package com.hartwig.hmftools.iclusion.api;

import com.squareup.moshi.Json;

class IclusionObjectGene {

    @Json(name = "id")
    public String id;

    @Json(name = "gene_name")
    public String geneName;
}
