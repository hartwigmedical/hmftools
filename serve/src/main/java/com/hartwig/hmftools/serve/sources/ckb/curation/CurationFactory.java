package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.Map;

import com.google.common.collect.Maps;

public class CurationFactory {

    static final Map<String, String> GENE_MAPPINGS = Maps.newHashMap();

    static final Map<CurationEntry, CurationEntry> MUTATION_MAPPINGS = Maps.newHashMap();

    private CurationFactory() {
    }

    static {

        MUTATION_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry("-", "MSI HIGH"));
        MUTATION_MAPPINGS.put(new CurationEntry("MSI", "low"), new CurationEntry("-", "MSI LOW"));
        MUTATION_MAPPINGS.put(new CurationEntry("TMB", "high"), new CurationEntry("-", "TumMutLoad HIGH"));
        MUTATION_MAPPINGS.put(new CurationEntry("TMB", "low"), new CurationEntry("-", "TumMutLoad LOW"));

    }
}
