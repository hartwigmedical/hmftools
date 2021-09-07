package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

final class DrugCurationFactory {

    static final Map<DrugCurationKey, DrugCurationValues> DRUG_MAPPINGS = Maps.newHashMap();

    static final Set<DrugCurationKey> DRUG_BLACKLIST = Sets.newHashSet();

    private DrugCurationFactory() {
    }

    static {
        populateBlacklist();

        populateRenames();
        populateAndMappings();
        populateOrMappings();
    }

    private static void populateBlacklist() {
        DRUG_BLACKLIST.add(civicB("119"));

        // The below drugs are supposed to be given sequentially.
        DRUG_BLACKLIST.add(civicB("Sunitinib,Everolimus"));
        DRUG_BLACKLIST.add(civicB("Vismodegib,Erismodegib"));
    }

    private static void populateRenames() {
        DRUG_MAPPINGS.put(cgiA("Cetuximab,Mab"), rename("Cetuximab"));
        DRUG_MAPPINGS.put(cgiA("Egfr"), rename("EGFR inhibitor"));
        DRUG_MAPPINGS.put(cgiA("Panitumumab,Mab"), rename("Panitumumab"));
        DRUG_MAPPINGS.put(cgiA("Pertuzumab,Mab"), rename("Pertuzumab"));
        DRUG_MAPPINGS.put(cgiA("Trastuzumab,Mab"), rename("Trastuzumab"));

        DRUG_MAPPINGS.put(cgiB("1st"), rename("EGFR inhibitor (1st gen)"));
        DRUG_MAPPINGS.put(cgiB("2'-deoxyinosine"), rename("EGFR inhibitor (2nd gen)"));
        DRUG_MAPPINGS.put(cgiB("2-butanone"), rename("MEK inhibitor"));
        DRUG_MAPPINGS.put(cgiB("Cetuximab,Mab"), rename("Cetuximab"));
        DRUG_MAPPINGS.put(cgiB("Egfr"), rename("EGFR inhibitor"));
        DRUG_MAPPINGS.put(cgiB("Mab"), rename("Anti-EGFR monoclonal antibody"));
        DRUG_MAPPINGS.put(cgiB("Panitumumab,Mab"), rename("Panitumumab"));
        DRUG_MAPPINGS.put(cgiB("Platinum"), rename("Platinum agent"));

        DRUG_MAPPINGS.put(civicB("120"), rename("AG-120"));
        DRUG_MAPPINGS.put(civicB("398"), rename("BGJ398"));
        DRUG_MAPPINGS.put(civicB("Alpelisib,Nvp"), rename("Alpelisib"));
        DRUG_MAPPINGS.put(civicB("Alpha"), rename("Pegylated IFN-alpha-2a"));
        DRUG_MAPPINGS.put(civicB("Anti"), rename("Anti-EGFR monoclonal antibody"));
        DRUG_MAPPINGS.put(civicB("Azd,5363"), rename("AZD5363"));
        DRUG_MAPPINGS.put(civicB("Azd-9291"), rename("Osimertinib"));
        DRUG_MAPPINGS.put(civicB("Azd4547"), rename("AZD4547"));
        DRUG_MAPPINGS.put(civicB("Ci-1040"), rename("CI-1040"));
        DRUG_MAPPINGS.put(civicB("Egfr"), rename("EGFR inhibitor"));
        DRUG_MAPPINGS.put(civicB("Folfoxprotocol"), rename("FOLFOX regimen"));
        DRUG_MAPPINGS.put(civicB("Jnj"), rename("Erdafitinib"));
        DRUG_MAPPINGS.put(civicB("Loxo,101"), rename("Larotrectinib"));
        DRUG_MAPPINGS.put(civicB("Mtor inhibitors"), rename("mTOR inhibitors"));
        DRUG_MAPPINGS.put(civicB("Nvp-bgj398"), rename("BGJ398"));
        DRUG_MAPPINGS.put(civicB("Nvp,Dactolisib"), rename("Dactolisib"));
        DRUG_MAPPINGS.put(civicB("Ro4987655"), rename("RO4987655"));
        DRUG_MAPPINGS.put(civicB("Rg7112"), rename("RG7112"));
        DRUG_MAPPINGS.put(civicB("Trans"), rename("Tretinoin"));
        DRUG_MAPPINGS.put(civicB("Vegf"), rename("Anti-VEGF monoclonal antibody"));
    }

    private static void populateAndMappings() {
        DRUG_MAPPINGS.put(cgiA("Dabrafenib,Trametinib"), and("Dabrafenib", "Trametinib"));
        DRUG_MAPPINGS.put(cgiA("Tretinoin,Arsenic"), and("Tretinoin", "Arsenic Trioxide"));
        DRUG_MAPPINGS.put(cgiA("Vemurafenib,Cobimetinib"), and("Vemurafenib", "Cobimetinib"));

        DRUG_MAPPINGS.put(cgiB("Everolimus,Trastuzumab,Mab"), and("Everolimus", "Trastuzumab", "Chemotherapy"));
        DRUG_MAPPINGS.put(cgiB("Trastuzumab,Lapatinib,Mab"), and("ETrastuzumab", "Lapatinib"));

        DRUG_MAPPINGS.put(civicA("Dabrafenib,Trametinib"), and("Dabrafenib", "Trametinib"));
        DRUG_MAPPINGS.put(civicA("Docetaxel,Selumetinib"), and("Docetaxel", "Selumetinib"));
        DRUG_MAPPINGS.put(civicA("Trametinib,Dabrafenib"), and("Trametinib", "Dabrafenib"));
        DRUG_MAPPINGS.put(civicA("Trans,Arsenic"), and("Tretinoin", "Arsenic Trioxide"));

        DRUG_MAPPINGS.put(civicB("Afatinib,Trastuzumab"), and("Afatinib", "Trastuzumab"));
        DRUG_MAPPINGS.put(civicB("Bevacizumab,Erlotinib"), and("Bevacizumab", "Erlotinib"));
        DRUG_MAPPINGS.put(civicB("Binimetinib,Ribociclib"), and("Binimetinib", "Ribociclib"));
        DRUG_MAPPINGS.put(civicB("Buparlisib,Carboplatin,Paclitaxel"), and("Buparlisib", "Carboplatin", "Paclitaxel"));
        DRUG_MAPPINGS.put(civicB("Carboplatin,Paclitaxel"), and("Carboplatin", "Paclitaxel"));
        DRUG_MAPPINGS.put(civicB("Cediranib,Olaparib"), and("Cediranib", "Olaparib"));
        DRUG_MAPPINGS.put(civicB("Cisplatin,Gemcitabine"), and("Cisplatin", "Gemcitabine"));
        DRUG_MAPPINGS.put(civicB("Cisplatin,Vinorelbine"), and("Cisplatin", "Vinorelbine"));
        DRUG_MAPPINGS.put(civicB("Cobimetinib,Vemurafenib"), and("Cobimetinib", "Vemurafenib"));
        DRUG_MAPPINGS.put(civicB("Cytarabine,Daunorubicin"), and("Cytarabine", "Daunorubicin"));
        DRUG_MAPPINGS.put(civicB("Dabrafenib,Trametinib"), and("Dabrafenib", "Trametinib"));
        DRUG_MAPPINGS.put(civicB("Docetaxel,Carboplatin,Sorafenib"), and("Docetaxel", "Carboplatin", "Sorafenib"));
        DRUG_MAPPINGS.put(civicB("Docetaxel,Cetuximab"), and("Docetaxel", "Cetuximab"));
        DRUG_MAPPINGS.put(civicB("Docetaxel,Selumetinib"), and("Docetaxel", "Selumetinib"));
        DRUG_MAPPINGS.put(civicB("Docetaxel,Trametinib"), and("Docetaxel", "Trametinib"));
        DRUG_MAPPINGS.put(civicB("Encorafenib,Cetuximab"), and("Encorafenib", "Cetuximab"));
        DRUG_MAPPINGS.put(civicB("Erlotinib,Gemcitabine"), and("Erlotinib", "Gemcitabine"));
        DRUG_MAPPINGS.put(civicB("Fluorouracil,Oxaliplatin"), and("Fluorouracil", "Oxaliplatin"));
        DRUG_MAPPINGS.put(civicB("Lapatinib,Capecitabine"), and("Lapatinib", "Capecitabine"));
        DRUG_MAPPINGS.put(civicB("Lapatinib,Trastuzumab"), and("Lapatinib", "Trastuzumab"));
        DRUG_MAPPINGS.put(civicB("Palbociclib,Letrozole"), and("Palbociclib", "Letrozole"));
        DRUG_MAPPINGS.put(civicB("Pertuzumab,Trastuzumab,Docetaxel"), and("Pertuzumab", "Trastuzumab", "Docetaxel"));
        DRUG_MAPPINGS.put(civicB("Refametinib,Sorafenib"), and("Refametinib", "Sorafenib"));
        DRUG_MAPPINGS.put(civicB("Selumetinib,Irinotecan"), and("Selumetinib", "Irinotecan"));
        DRUG_MAPPINGS.put(civicB("Sorafenib,Paclitaxel,Carboplatin"), and("Sorafenib", "Paclitaxel", "Carboplatin"));
        DRUG_MAPPINGS.put(civicB("Taxane,Platinum"), and("Taxane", "Platinum"));
        DRUG_MAPPINGS.put(civicB("Trametinib,Docetaxel,Pemetrexed"), and("Trametinib", "Docetaxel", "Pemetrexed"));
        DRUG_MAPPINGS.put(civicB("Trastuzumab,Capecitabine"), and("Trastuzumab", "Capecitabine"));
        DRUG_MAPPINGS.put(civicB("Trastuzumab,Lapatinib"), and("Trastuzumab", "Lapatinib"));
        DRUG_MAPPINGS.put(civicB("Vemurafenib,Cetuximab,Irinotecan"), and("Vemurafenib", "Cetuximab", "Irinotecan"));
        DRUG_MAPPINGS.put(civicB("Vemurafenib,Cobimetinib"), and("Vemurafenib", "Cobimetinib"));
    }

    private static void populateOrMappings() {
        DRUG_MAPPINGS.put(cgiA("Nilotinib,Dasatinib"), or("Nilotinib", "Dasatinib"));

        DRUG_MAPPINGS.put(civicA("Cetuximab,Panitumumab"), or("Cetuximab", "Panitumumab"));
        DRUG_MAPPINGS.put(civicA("Dasatinib,Nilotinib"), or("Dasatinib", "Nilotinib"));
        DRUG_MAPPINGS.put(civicA("Panitumumab,Cetuximab"), or("Panitumumab", "Cetuximab"));

        DRUG_MAPPINGS.put(civicB("Afatinib,Erlotinib,Gefitinib"), or("Afatinib", "Erlotinib", "Gefitinib"));
        DRUG_MAPPINGS.put(civicB("Afatinib,Lapatinib,Trastuzumab"), or("Afatinib", "Lapatinib", "Trastuzumab"));
        DRUG_MAPPINGS.put(civicB("Atezolizumab,Pembrolizumab,Nivolumab"), or("Atezolizumab", "Pembrolizumab", "Nivolumab"));
        DRUG_MAPPINGS.put(civicB("Carboplatin,Cisplatin"), or("Carboplatin", "Cisplatin"));
        DRUG_MAPPINGS.put(civicB("Cetuximab,Panitumumab"), or("Cetuximab", "Panitumumab"));
        DRUG_MAPPINGS.put(civicB("Cisplatin,Carboplatin"), or("Cisplatin", "Carboplatin"));
        DRUG_MAPPINGS.put(civicB("Dacomitinib,Erlotinib"), or("Dacomitinib", "Erlotinib"));
        DRUG_MAPPINGS.put(civicB("Docetaxel,Vinorelbine,Gemcitabine"), or("Docetaxel", "Vinorelbine", "Gemcitabine"));
        DRUG_MAPPINGS.put(civicB("Erlotinib,Gefitinib"), or("Erlotinib", "Gefitinib"));
        DRUG_MAPPINGS.put(civicB("Gefitinib,Erlotinib,Afatinib"), or("Gefitinib", "Erlotinib", "Afatinib"));
        DRUG_MAPPINGS.put(civicB("Gemcitabine,Vinorelbine,Paclitaxel"), or("Gemcitabine", "Vinorelbine", "Paclitaxel"));
        DRUG_MAPPINGS.put(civicB("Nilutamide,Cyproterone,Bornyl acetate,Flutamide,Bicalutamide"),
                or("Nilutamide", "Cyproterone", "Bornyl acetate", "Flutamide", "Bicalutamide"));
        DRUG_MAPPINGS.put(civicB("Nivolumab,Atezolizumab"), or("Nivolumab", "Atezolizumab"));
        DRUG_MAPPINGS.put(civicB("Panitumumab,Cetuximab"), or("Panitumumab", "Cetuximab"));
        DRUG_MAPPINGS.put(civicB("Quizartinib,Sorafenib"), or("Quizartinib", "Sorafenib"));
        DRUG_MAPPINGS.put(civicB("Vinorelbine,Docetaxel,Gemcitabine"), or("Vinorelbine", "Docetaxel", "Gemcitabine"));
    }

    @NotNull
    private static DrugCurationValues rename(@NotNull String drug) {
        List<List<String>> drugs = Lists.newArrayList();
        drugs.add(Lists.newArrayList(drug));
        return ImmutableDrugCurationValues.builder().drugs(drugs).build();
    }

    @NotNull
    private static DrugCurationValues and(@NotNull String... drugs) {
        if (drugs.length < 2) {
            throw new IllegalStateException("At least 2 drugs required for an AND relation!");
        }
        List<List<String>> curatedDrugs = Lists.newArrayList();
        curatedDrugs.add(Lists.newArrayList(drugs));
        return ImmutableDrugCurationValues.builder().drugs(curatedDrugs).build();
    }

    @NotNull
    private static DrugCurationValues or(@NotNull String... drugs) {
        if (drugs.length < 2) {
            throw new IllegalStateException("At least 2 drugs required for an OR relation!");
        }
        List<List<String>> curatedDrugs = Lists.newArrayList();
        for (String drug : drugs) {
            curatedDrugs.add(Lists.newArrayList(drug));
        }
        return ImmutableDrugCurationValues.builder().drugs(curatedDrugs).build();
    }

    @NotNull
    private static DrugCurationKey cgiA(@NotNull String treatment) {
        return new DrugCurationKey(ViccSource.CGI, EvidenceLevel.A, treatment);
    }

    @NotNull
    private static DrugCurationKey cgiB(@NotNull String treatment) {
        return new DrugCurationKey(ViccSource.CGI, EvidenceLevel.B, treatment);
    }

    @NotNull
    private static DrugCurationKey civicA(@NotNull String treatment) {
        return new DrugCurationKey(ViccSource.CIVIC, EvidenceLevel.A, treatment);
    }

    @NotNull
    private static DrugCurationKey civicB(@NotNull String treatment) {
        return new DrugCurationKey(ViccSource.CIVIC, EvidenceLevel.B, treatment);
    }
}
