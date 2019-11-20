package com.hartwig.hmftools.vicc.reader;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

final class ViccDatamodelCheckerFactory {

    private ViccDatamodelCheckerFactory() {
    }

    @NotNull
    static ViccDatamodelChecker viccEntryDatamodelChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("source", true);
        map.put("genes", true);
        map.put("gene_identifiers", true);
        map.put("feature_names", false);
        map.put("features", true);
        map.put("association", true);
        map.put("tags", true);
        map.put("dev_tags", true);
        map.put("cgi", false);
        map.put("brca", false);
        map.put("sage", false);
        map.put("pmkb", false);
        map.put("oncokb", false);
        map.put("jax", false);
        map.put("jax_trials", false);
        map.put("molecularmatch", false);
        map.put("molecularmatch_trials", false);
        map.put("civic", false);
        return new ViccDatamodelChecker("ViccEntry", map);
    }
}
