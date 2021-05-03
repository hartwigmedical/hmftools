package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.Map;

import com.google.common.collect.Maps;

public class CurationFactory {

    static final Map<CurationEntry, CurationEntry> VARIANT_MAPPINGS = Maps.newHashMap();

    static final Map<String, String> GENE_MAPPINGS = Maps.newHashMap();

    private CurationFactory() {
    }

    static {
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry("-", "MSI high"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "low"), new CurationEntry("-", "MSI low"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "negative"), new CurationEntry("-", "MSI neg"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "high"), new CurationEntry("-", "TMB high"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "low"), new CurationEntry("-", "TMB low"));

        GENE_MAPPINGS.put("SEPTIN14", "SEPT14");
    }
}
