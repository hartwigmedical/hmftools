package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import com.squareup.moshi.Json;

class IclusionObjectStudy {
    @Json(name = "id") public String id;
    @Json(name = "title") public String title;
    @Json(name = "acronym") public String acronym;
    @Json(name = "eudra") public String eudra;
    @Json(name = "nct") public String nct;
    @Json(name = "ipn") public String ipn;
    @Json(name = "ccmo") public String ccmo;
    @Json(name = "type_name")  public String typeName;
    @Json(name = "age_name")  public String ageName;
    @Json(name = "phase_name")  public String phaseName;
    @Json(name = "indication_ids") public List<String> indicationIds;
    @Json(name = "mutations") public List<IclusionObjectMutation> mutations;
}
