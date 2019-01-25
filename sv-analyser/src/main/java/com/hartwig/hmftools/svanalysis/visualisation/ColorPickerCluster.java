package com.hartwig.hmftools.svanalysis.visualisation;

import static java.util.stream.Collectors.toMap;

import java.awt.Color;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class ColorPickerCluster implements ColorPicker {


    private static final Color[] COLOURS = new Color[] {COLOR2, COLOR3, COLOR4, COLOR5, COLOR6, COLOR7, COLOR8, COLOR9, COLOR10, COLOR11 };

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

            String color = i < COLOURS.length ? toString(COLOURS[i]) : "black";
            if (clusterSize.count == 1) {
                color = toString(COLOR1);
            }

            clusterIdMap.put(clusterSizeList.get(i).clusterId, color);
        }

    }

    public static String toString(Color color) {
        return color.getRed() + "," + color.getGreen() + "," + color.getBlue();
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
