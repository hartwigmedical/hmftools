package com.hartwig.hmftools.serve.sources.docm.curation;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

final class CurationFactory {

    static final Set<CurationKey> ENTRY_BLACKLIST = Sets.newHashSet();

    static final Map<String, String> GENE_MAPPINGS = Maps.newHashMap();

    static {
        // Some gene names in DoCM are not HGNC-compliant, so we map them to the HGNC versions
        GENE_MAPPINGS.put("H3F3A", "H3-3A");
        GENE_MAPPINGS.put("RQCD1", "CNOT9");
        GENE_MAPPINGS.put("TCEB1", "ELOC");
        GENE_MAPPINGS.put("HIST1H3I", "K28M");

        // Not clear what the "minus" means, so ignoring. Could be DEL?
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000288602", "K601-"));
        ENTRY_BLACKLIST.add(new CurationKey("CTNNB1", "ENST00000349496", "VSHWQQQSYLDSGIHSG22-"));
        ENTRY_BLACKLIST.add(new CurationKey("CTNNB1", "ENST00000349496", "WQQQSYLD25-"));
        ENTRY_BLACKLIST.add(new CurationKey("FLT3", "ENST00000380982", "D835-"));
        ENTRY_BLACKLIST.add(new CurationKey("FLT3", "ENST00000380982", "I836-"));
        ENTRY_BLACKLIST.add(new CurationKey("KIT", "ENST00000288135", "V559-"));
        ENTRY_BLACKLIST.add(new CurationKey("PDGFRA", "ENST00000257290", "RD841-"));
        ENTRY_BLACKLIST.add(new CurationKey("PDGFRA", "ENST00000257290", "DIM842-"));
        ENTRY_BLACKLIST.add(new CurationKey("PDGFRA", "ENST00000257290", "DIMH842-"));
        ENTRY_BLACKLIST.add(new CurationKey("PDGFRA", "ENST00000257290", "IMHD843-"));
        ENTRY_BLACKLIST.add(new CurationKey("PIK3R1", "ENST00000396611", "T576-"));
        ENTRY_BLACKLIST.add(new CurationKey("PIK3R1", "ENST00000396611", "DKRMNS560-"));
        ENTRY_BLACKLIST.add(new CurationKey("PTEN", "ENST00000371953", "K267X"));
        ENTRY_BLACKLIST.add(new CurationKey("RET", "ENST00000355710", "FPEEEKCFC612-"));
        ENTRY_BLACKLIST.add(new CurationKey("RET", "ENST00000355710", "EL632-"));
        ENTRY_BLACKLIST.add(new CurationKey("RET", "ENST00000355710", "DVYE898-"));

        // Unclear what these variants means. Maybe these are INS?
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000288602", "GD593GN"));
        ENTRY_BLACKLIST.add(new CurationKey("ERBB2", "ENST00000269571", "G776CX"));
        ENTRY_BLACKLIST.add(new CurationKey("ERBB2", "ENST00000269571", "G778GSPX"));
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288SX"));
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288PX"));
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288HX"));
        ENTRY_BLACKLIST.add(new CurationKey("WT1", "ENST00000332351", "S381LYG"));

        // Wildcard variants are ignored
        ENTRY_BLACKLIST.add(new CurationKey("APC", "ENST00000257430", "KIGT1310X"));
        ENTRY_BLACKLIST.add(new CurationKey("APC", "ENST00000257430", "KIGTRSA1310X"));
        ENTRY_BLACKLIST.add(new CurationKey("APC", "ENST00000257430", "GP1466X"));
        ENTRY_BLACKLIST.add(new CurationKey("APC", "ENST00000257430", "IDS1557X"));
        ENTRY_BLACKLIST.add(new CurationKey("ATR", "ENST00000350721", "I774X"));
        ENTRY_BLACKLIST.add(new CurationKey("ERBB2", "ENST00000269571", "L755X"));
        ENTRY_BLACKLIST.add(new CurationKey("HERC2", "ENST00000261609", "D759X"));
        ENTRY_BLACKLIST.add(new CurationKey("KIT", "ENST00000288135", "V560X"));
        ENTRY_BLACKLIST.add(new CurationKey("KIT", "ENST00000288135", "P577X"));
        ENTRY_BLACKLIST.add(new CurationKey("PIK3R1", "ENST00000396611", "E439X"));
        ENTRY_BLACKLIST.add(new CurationKey("PIK3R1", "ENST00000396611", "W583X"));

        // Not sure what these mean? Maybe wildcards?
        ENTRY_BLACKLIST.add(new CurationKey("ABCB1", "ENST00000265724", "I1145"));
        ENTRY_BLACKLIST.add(new CurationKey("ETS2", "ENST00000360214", "P341"));
        ENTRY_BLACKLIST.add(new CurationKey("KIT", "ENST00000288135", "L862"));
        ENTRY_BLACKLIST.add(new CurationKey("MGMT", "ENST00000306010", "R22"));

        // Variants that don't exist on the configured transcript
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "F203L"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "T207I"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "G77L"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "G77V"));
    }

    private CurationFactory() {
    }
}
