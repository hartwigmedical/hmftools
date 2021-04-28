package com.hartwig.hmftools.common.genome.refgenome;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class GeneNameMapping {

    @NotNull
    private final Map<String, String> v37Map;
    @NotNull
    private final Map<String, String> v38Map;

    @NotNull
    public static GeneNameMapping loadFromEmbeddedResource() {
        final Map<String, String> v37Map = Maps.newHashMap();
        final Map<String, String> v38Map = Maps.newHashMap();
        final InputStream inputStream = GeneNameMapping.class.getResourceAsStream("/refgenome/gene_name_mapping.tsv");

        // Skip header
        new BufferedReader(new InputStreamReader(inputStream)).lines().skip(1)
                .map(x -> x.split("\t"))
                .forEach(strings -> {
                    String v37 = strings[0];
                    String v38 = strings[1];
                    if (!v38.equals("NA")) {
                        v38Map.put(v38, v37);
                    }
                    v37Map.put(v37, v38);
                });

        return new GeneNameMapping(v37Map, v38Map);
    }

    public GeneNameMapping(@NotNull final Map<String, String> v37Map, @NotNull final Map<String, String> v38Map) {
        this.v37Map = v37Map;
        this.v38Map = v38Map;
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
    public String v38Gene(@NotNull final String v37GeneId) {
        String result = v37Map.get(v37GeneId);
        if (result == null ) {
            return "NA";
        }
        return result;
    }
}
