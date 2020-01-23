package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.squareup.moshi.Json;

public class IclusionGene {
    @Json(name = "id") public String id;
    @Json(name = "gene_name") public String geneName;
    @Json(name = "variant_ids") public List<String> variantIds;
}
