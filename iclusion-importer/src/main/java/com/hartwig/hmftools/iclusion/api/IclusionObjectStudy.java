package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.squareup.moshi.Json;

class IclusionObjectStudy {
    @Json(name = "id")
    public String id;
    @Json(name = "title")
    public String title;
    @Json(name = "acronym")
    public String acronym;
    @Json(name = "eudra")
    public String eudra;
    @Json(name = "nct")
    public String nct;
    @Json(name = "ipn")
    public String ipn;
    @Json(name = "ccmo")
    public String ccmo;
    @Json(name = "indication_ids")
    public List<String> indicationIds;
    @Json(name = "blacklist_indication_ids")
    public List<String> blacklistIndicationIds;
    @Json(name = "mutation_conditions")
    public List<IclusionObjectMutationCondition> mutationConditions;
}
