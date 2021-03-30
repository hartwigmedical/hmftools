package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class DrugCurationFactory {

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
    }

    private static void populateRenames() {

    }

    private static void populateAndMappings() {

    }

    private static void populateOrMappings() {
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

}
