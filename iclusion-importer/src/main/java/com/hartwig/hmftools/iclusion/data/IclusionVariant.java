package com.hartwig.hmftools.iclusion.data;

import java.util.List;

import com.squareup.moshi.Json;

public class IclusionVariant {
    @Json(name = "id") public String id;
    @Json(name = "variant_name") public String variantName;
    @Json(name = "gene_ids") public List<String> geneIds;
}
