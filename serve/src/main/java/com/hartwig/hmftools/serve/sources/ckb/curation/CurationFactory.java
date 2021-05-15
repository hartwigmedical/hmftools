package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.ckb.classification.CkbConstants;

public class CurationFactory {

    static final Map<CurationEntry, CurationEntry> VARIANT_MAPPINGS = Maps.newHashMap();

    private CurationFactory() {
    }

    static {
        // CKB uses "genes" to model evidence on characteristics. We map this away from genes.
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry(CkbConstants.NO_GENE, "MSI high"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "low"), new CurationEntry(CkbConstants.NO_GENE, "MSI low"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "negative"), new CurationEntry(CkbConstants.NO_GENE, "MSI neg"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "high"), new CurationEntry(CkbConstants.NO_GENE, "TMB high"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "low"), new CurationEntry(CkbConstants.NO_GENE, "TMB low"));
    }
}
