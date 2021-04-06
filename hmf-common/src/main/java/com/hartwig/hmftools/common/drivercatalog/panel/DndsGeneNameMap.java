package com.hartwig.hmftools.common.drivercatalog.panel;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class DndsGeneNameMap {

    private final Map<String, String> v38Map = Maps.newHashMap();
    private final Map<String, String> v37Map = Maps.newHashMap();

    public DndsGeneNameMap() {
        final InputStream inputStream = DriverGenePanel.class.getResourceAsStream("/drivercatalog/DndsGeneNameMap.tsv");
        new BufferedReader(new InputStreamReader(inputStream)).lines()
                .filter(x -> !x.startsWith("v37"))
                .map(x -> x.split("\t"))
                .forEach(strings -> {
                    String v37 = strings[0];
                    String v38 = strings[1];
                    if (!v38.equals("NA")) {
                        v38Map.put(v38, v37);
                    }
                    v37Map.put(v37, v38);
                });
    }

    public boolean isValidV38Gene(@NotNull final String v38GeneId) {
        return v38Map.containsKey(v38GeneId) && !v38Map.get(v38GeneId).equals("NA");
    }

    public boolean isValidV37Gene(@NotNull final String v37GeneId) {
        return v37Map.containsKey(v37GeneId);
    }

    @NotNull
    public String v37Gene(@NotNull final String v38GeneId) {
        String result = v38Map.get(v38GeneId);
        if (result == null || result.equals("NA")) {
            throw new IllegalArgumentException("Invalid v38 gene " + v38GeneId);
        }
        return result;
    }

    @NotNull
    String v38Gene(@NotNull final String v37GeneId) {
        String result = v37Map.get(v37GeneId);
        if (result == null ) {
            return "NA";
        }
        return result;
    }
}
