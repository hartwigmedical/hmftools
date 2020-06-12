package com.hartwig.hmftools.serve.vicc.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

final class CurationFactory {

    static final Map<CurationKey, String> ONCOKB_FEATURE_NAME_MAPPINGS = Maps.newHashMap();

    static final Set<CurationKey> ONCOKB_FEATURE_BLACKLIST = Sets.newHashSet();

    static {
        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("EPAS1", "ENST00000263734", "533_534del"), "I533_P534del");
        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("KIT", "ENST00000288135", "V559del"),"V560del");
        ONCOKB_FEATURE_NAME_MAPPINGS.put(new CurationKey("PTEN", "ENST00000371953", "I32del"), "I33del");

        //        ONCOKB_FEATURE_BLACKLIST.add(new FilterKey("RAD"))
    }

    private CurationFactory() {
    }

}
