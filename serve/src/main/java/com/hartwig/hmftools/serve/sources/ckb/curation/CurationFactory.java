package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.Map;

import com.google.common.collect.Maps;

public class CurationFactory {

    static final Map<CurationEntry, CurationEntry> VARIANT_MAPPINGS = Maps.newHashMap();

    static final Map<String, String> GENE_MAPPINGS = Maps.newHashMap();

    private CurationFactory() {
    }

    static {
        // CKB uses "genes" to model evidence on characteristics. We map this away from genes.
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "high"), new CurationEntry("-", "MSI high"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "low"), new CurationEntry("-", "MSI low"));
        VARIANT_MAPPINGS.put(new CurationEntry("MSI", "negative"), new CurationEntry("-", "MSI neg"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "high"), new CurationEntry("-", "TMB high"));
        VARIANT_MAPPINGS.put(new CurationEntry("TMB", "low"), new CurationEntry("-", "TMB low"));

        // These are all valid gene symbols, for which we use a different synonym compared to CKB.
        GENE_MAPPINGS.put("ABRAXAS1", "FAM175A");
        GENE_MAPPINGS.put("ACP3", "ACPP");
        GENE_MAPPINGS.put("CARS1", "CARS");
        GENE_MAPPINGS.put("CEP43", "FGFR1OP");
        GENE_MAPPINGS.put("CIP2A", "KIAA1524");
        GENE_MAPPINGS.put("CCN2", "CTGF");
        GENE_MAPPINGS.put("CCN6", "WISP3");
        GENE_MAPPINGS.put("COP1", "RFWD2");
        GENE_MAPPINGS.put("FYB1", "FYB");
        GENE_MAPPINGS.put("GSDME", "DFNA5");
        GENE_MAPPINGS.put("H1-2", "HIST1H1C");
        GENE_MAPPINGS.put("H3-3A", "H3F3A");
        GENE_MAPPINGS.put("H3-3B", "H3F3B");
        GENE_MAPPINGS.put("H3-4", "HIST3H3");
        GENE_MAPPINGS.put("H3-5", "H3F3C");
        GENE_MAPPINGS.put("H2BC4", "HIST1H2BC");
        GENE_MAPPINGS.put("H2BC5", "HIST1H2BC");
        GENE_MAPPINGS.put("H3C1", "HIST1H3A");
        GENE_MAPPINGS.put("H3C13", "HIST2H3D");
        GENE_MAPPINGS.put("H3C2", "HIST1H3B");
        GENE_MAPPINGS.put("H3C3", "HIST1H3C");
        GENE_MAPPINGS.put("H3C4", "HIST1H3D");
        GENE_MAPPINGS.put("H3C6", "HIST1H3E");
        GENE_MAPPINGS.put("H3C7", "HIST1H3F");
        GENE_MAPPINGS.put("H3C8", "HIST1H3G");
        GENE_MAPPINGS.put("H3C10", "HIST1H3H");
        GENE_MAPPINGS.put("H3C11", "HIST1H3I");
        GENE_MAPPINGS.put("H3C12", "HIST1H3J");
        GENE_MAPPINGS.put("H3C14", "HIST2H3C");
        GENE_MAPPINGS.put("H3C15", "HIST2H3A");
        GENE_MAPPINGS.put("OGA", "MGEA5");
        GENE_MAPPINGS.put("RARS1", "RARS");
        GENE_MAPPINGS.put("RELCH", "KIAA1468");
        GENE_MAPPINGS.put("SEPTIN3", "SEPT3");
        GENE_MAPPINGS.put("SEPTIN6", "SEPT6");
        GENE_MAPPINGS.put("SEPTIN7", "SEPT7");
        GENE_MAPPINGS.put("SEPTIN9", "SEPT9");
        GENE_MAPPINGS.put("SEPTIN14", "SEPT14");
        GENE_MAPPINGS.put("SLC66A1", "PQLC2");
        GENE_MAPPINGS.put("STING1", "TMEM173");
        GENE_MAPPINGS.put("TENT5C", "FAM46C");
        GENE_MAPPINGS.put("TBXT", "T");
        GENE_MAPPINGS.put("ZFTA", "C11orf95");
    }
}
