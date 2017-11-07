package com.hartwig.hmftools.common.centromeres;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

public class Centromeres {

    public static Map<String, GenomeRegion> grch37() {
        final Map<String, GenomeRegion> result = Maps.newHashMap();
        addCentromere(result, "1", 121535434, 124535434);
        addCentromere(result, "2", 92326171, 95326171);
        addCentromere(result, "3", 90504854, 93504854);
        addCentromere(result, "4", 49660117, 52660117);
        addCentromere(result, "5", 46405641, 49405641);
        addCentromere(result, "6", 58830166, 61830166);
        addCentromere(result, "7", 58054331, 61054331);
        addCentromere(result, "8", 43838887, 46838887);
        addCentromere(result, "9", 47367679, 50367679);
        addCentromere(result, "10", 39254935, 42254935);
        addCentromere(result, "11", 51644205, 54644205);
        addCentromere(result, "12", 34856694, 37856694);
        addCentromere(result, "13", 16000000, 19000000);
        addCentromere(result, "14", 16000000, 19000000);
        addCentromere(result, "15", 17000000, 20000000);
        addCentromere(result, "16", 35335801, 38335801);
        addCentromere(result, "17", 22263006, 25263006);
        addCentromere(result, "18", 15460898, 18460898);
        addCentromere(result, "19", 24681782, 27681782);
        addCentromere(result, "20", 26369569, 29369569);
        addCentromere(result, "21", 11288129, 14288129);
        addCentromere(result, "22", 13000000, 16000000);
        addCentromere(result, "X", 58632012, 61632012);
        addCentromere(result, "Y", 10104553, 13104553);

        return result;
    }

    private static void addCentromere(@NotNull final Map<String, GenomeRegion> map, @NotNull final String chromsome, final long start,
            final long end) {
        map.put(chromsome, GenomeRegionFactory.create(chromsome, start, end));
    }
}
