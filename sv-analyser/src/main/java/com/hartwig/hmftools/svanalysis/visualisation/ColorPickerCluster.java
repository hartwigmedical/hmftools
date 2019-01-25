package com.hartwig.hmftools.svanalysis.visualisation;

import static java.util.stream.Collectors.toMap;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class ColorPickerCluster implements ColorPicker {

    private static final String[] COLOURS = new String[] {
            //"166,206,227",
            //"31,120,180",
            "51,160,44", "227,26,28", "255,127,0", "202,178,214", "106,61,154", "178,223,138", "251,154,153", "253,191,111",
            "255,255,153" };

    private final Map<Integer, String> clusterIdMap = Maps.newHashMap();

    public ColorPickerCluster(@NotNull final List<Link> links) {
        final Comparator<ClusterSize> longComparator = Comparator.<ClusterSize>comparingLong(x -> x.count).reversed();

        final List<ClusterSize> clusterSizeList = Lists.newArrayList();
        links.stream()
                .map(Link::clusterId)
                .collect(toMap(x -> x, (x) -> 1L, Math::addExact))
                .forEach((key, value) -> clusterSizeList.add(new ClusterSize(key, value)));

        clusterSizeList.sort(longComparator);

        for (int i = 0; i < clusterSizeList.size(); i++) {
            final ClusterSize clusterSize = clusterSizeList.get(i);

            String color = i < COLOURS.length ? COLOURS[i] : "black";
            if (clusterSize.count == 1) {
                color = "31,120,180";
            }

            clusterIdMap.put(clusterSizeList.get(i).clusterId, color);
        }

    }

    @NotNull
    @Override
    public String color(final int clusterId, final int chainId) {
        return "color=" + clusterIdMap.get(clusterId);
    }

    private class ClusterSize {

        final int clusterId;
        final long count;

        ClusterSize(final int clusterId, final long count) {
            this.clusterId = clusterId;
            this.count = count;
        }
    }

}
