package com.hartwig.hmftools.common.chromosome;

import java.util.Map;

import com.google.common.collect.Maps;

@Deprecated
public final class ChromosomeLengths {

    private ChromosomeLengths() {
    }

    public static Map<String, Integer> getChromosomeLengths() {
        Map<String, Integer> lengths = Maps.newHashMap();

        // taken from GRCH37
        lengths.put("1", 249250621);
        lengths.put("2", 243199373);
        lengths.put("3", 198022430);
        lengths.put("4", 191154276);
        lengths.put("5", 180915260);
        lengths.put("6", 171115067);
        lengths.put("7", 159138663);
        lengths.put("8", 146364022);
        lengths.put("9", 141213431);
        lengths.put("10", 135534747);
        lengths.put("11", 135006516);
        lengths.put("12", 133851895);
        lengths.put("13", 115169878);
        lengths.put("14", 107349540);
        lengths.put("15", 102531392);
        lengths.put("16", 90354753);
        lengths.put("17", 83257441);
        lengths.put("18", 78077248);
        lengths.put("19", 59128983);
        lengths.put("20", 63025520);
        lengths.put("21", 48129895);
        lengths.put("22", 51304566);
        lengths.put("X", 155270560);
        lengths.put("Y", 59373566);

        return lengths;
    }
}
