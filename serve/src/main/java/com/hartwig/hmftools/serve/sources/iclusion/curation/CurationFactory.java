package com.hartwig.hmftools.serve.sources.iclusion.curation;

import java.util.Map;

import com.google.common.collect.Maps;

final class CurationFactory {

    static final Map<String, String> GENE_MAPPINGS = Maps.newHashMap();

    static final Map<CurationEntry, CurationEntry> MUTATION_MAPPINGS = Maps.newHashMap();

    private CurationFactory() {
    }

    static {
        GENE_MAPPINGS.put("ERBB2 (HER2)", "ERBB2");
        GENE_MAPPINGS.put("PDGFRα", "PDGFRA");
        GENE_MAPPINGS.put("PDGFRβ", "PDGFRB");
        GENE_MAPPINGS.put("MRE11A", "MRE11"); //TODO: Removed when iClusion used correct gene name

        MUTATION_MAPPINGS.put(new CurationEntry("FGFR3", "FGFR3-WHSC1"), new CurationEntry("FGFR3", "FGFR3-NSD2"));
        MUTATION_MAPPINGS.put(new CurationEntry("MLL", "MLL-AF10"), new CurationEntry("KMT2A", "KMT2A-MLLT10"));
        MUTATION_MAPPINGS.put(new CurationEntry("MLL", "MLL-AF9"), new CurationEntry("KMT2A", "KMT2A-MLLT3"));

        MUTATION_MAPPINGS.put(new CurationEntry("MSI", "HIGH"), new CurationEntry("-", "MSI_HIGH"));
        MUTATION_MAPPINGS.put(new CurationEntry("TumMutLoad", "HIGH"), new CurationEntry("-", "TumMutLoad_HIGH"));
        MUTATION_MAPPINGS.put(new CurationEntry("HRD", "POSITIVE"), new CurationEntry("-", "HRD_POSITIVE"));
        MUTATION_MAPPINGS.put(new CurationEntry("HPV", "POSITIVE"), new CurationEntry("-", "HPV_POSITIVE"));
        MUTATION_MAPPINGS.put(new CurationEntry("EBV", "POSITIVE"), new CurationEntry("-", "EBV_POSITIVE"));
    }
}
