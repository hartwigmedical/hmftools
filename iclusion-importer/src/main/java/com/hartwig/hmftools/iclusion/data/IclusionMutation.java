package com.hartwig.hmftools.iclusion.data;

import com.squareup.moshi.Json;

public class IclusionMutation {
    @Json(name = "gene_id") public String geneId;
    @Json(name = "variant_id") public String variantId;
}
