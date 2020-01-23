package com.hartwig.hmftools.iclusion.data;

import java.util.List;

public class IclusionStudy {

    public String id;
    public String title;
    public String acronym;
    public String eudra;
    public String nct;
    public String ipn;
    public String ccmo;
    public List<String> indication_ids;
    public List<IclusionMutation> mutations;
}
