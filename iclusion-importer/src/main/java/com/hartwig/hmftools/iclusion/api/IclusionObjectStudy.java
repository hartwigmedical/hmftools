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
    @Json(name = "indication_ids") public List<String> indicationIds;
    @Json(name = "mutations") public List<IclusionObjectMutation> mutations;
    @Json(name = "type_id")  public String type_id;
    @Json(name = "type_name")  public String type_name;
    @Json(name = "age_id")  public String age_id;
    @Json(name = "age_name")  public String age_name;
    @Json(name = "phase_id")  public String phase_id;
    @Json(name = "phase_name")  public String phase_name;
}
