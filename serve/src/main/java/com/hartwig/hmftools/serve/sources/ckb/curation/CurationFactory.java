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
        VARIANT_MAPPINGS.put(new CurationEntry("HRD", "positive"), new CurationEntry(CkbConstants.NO_GENE, "HRD_pos"));
        VARIANT_MAPPINGS.put(new CurationEntry("HRD", "negative"), new CurationEntry(CkbConstants.NO_GENE, "HRD_neg"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry(CkbConstants.NO_GENE, "MSI_high"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry(CkbConstants.NO_GENE, "MSI_high"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "low"), new CurationEntry(CkbConstants.NO_GENE, "MSI_low"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "negative"), new CurationEntry(CkbConstants.NO_GENE, "MSI_neg"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "high"), new CurationEntry(CkbConstants.NO_GENE, "TMB_high"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "low"), new CurationEntry(CkbConstants.NO_GENE, "TMB_low"));
    }
}
