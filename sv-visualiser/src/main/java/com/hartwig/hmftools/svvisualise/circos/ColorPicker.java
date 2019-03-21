package com.hartwig.hmftools.svvisualise.circos;

import static java.util.stream.Collectors.toMap;

import java.awt.Color;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.svvisualise.data.Link;

import org.jetbrains.annotations.NotNull;

public class ColorPicker {

    private static final Color BLACK = new Color(1, 1, 1);

    private static final Color COLOR1 = new Color(31, 120, 180);
    private static final Color COLOR2 = new Color(227, 26, 28);
    private static final Color COLOR3 = new Color(255, 127, 0);
    private static final Color COLOR4 = new Color(51, 160, 44);
    private static final Color COLOR5 = new Color(106, 61, 154);
    private static final Color COLOR6 = new Color(140, 81, 10);
    private static final Color COLOR7 = new Color(1, 102, 94);
    private static final Color COLOR8 = new Color(212, 193, 23);

    private static final Color DEL = new Color(251, 154, 153);
    private static final Color DUP = new Color(178, 223, 138);
    private static final Color INS = new Color(255, 255, 153);
    private static final Color LINE = new Color(166,206,227);

    private static final Color[] COLOURS = new Color[] { COLOR1, COLOR2, COLOR3, COLOR4, COLOR5, COLOR6, COLOR7, COLOR8 };

    private final boolean clusterMode;
    private final Map<Integer, String> colorMap;
    private final double connectorTransparency;

    public ColorPicker(@NotNull final List<Link> links) {
        long clusterCount = links.stream().mapToLong(Link::clusterId).distinct().count();
        if (clusterCount > 1) {
            clusterMode = true;
            colorMap = clusterMap(links);
        } else {
            clusterMode = false;
            colorMap = chainMap(links);
        }

        connectorTransparency = connectorTransparency(links.size());
    }

    @NotNull
    public String connectorColor(final int clusterId, final int chainId) {
        String opaqueColor = color(clusterId, chainId);
        return opaqueColor.replace(")", "," + connectorTransparency + ")");
    }

    @NotNull
    public String color(final int clusterId, final int chainId) {
        if (clusterId == -1 || chainId == -1) {
            return toString(BLACK);
        }

        return clusterMode ? colorMap.get(clusterId) : colorMap.get(chainId);
    }

    @NotNull
    private static String toString(@NotNull final Color color) {
        return "color=(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ")";
    }

    @NotNull
    private static String simpleSvColor(@NotNull final String type) {
        switch (type) {
            case "DEL":
                return toString(DEL);
            case "DUP":
                return toString(DUP);
        }

        return toString(INS);
    }

    private class ClusterSize {

        final int clusterId;
        final long count;

        ClusterSize(final int clusterId, final long count) {
            this.clusterId = clusterId;
            this.count = count;
        }
    }

    @NotNull
    private Map<Integer, String> clusterMap(@NotNull final List<Link> links) {
        final Map<Integer, String> result = Maps.newHashMap();

        final Comparator<ClusterSize> longComparator = Comparator.<ClusterSize>comparingLong(x -> x.count).reversed();

        final List<ClusterSize> clusterSizeList = Lists.newArrayList();
        links.stream()
                .map(Link::clusterId)
                .collect(toMap(x -> x, (x) -> 1L, Math::addExact))
                .forEach((key, value) -> clusterSizeList.add(new ClusterSize(key, value)));

        clusterSizeList.sort(longComparator);

        for (int i = 0; i < clusterSizeList.size(); i++) {
            String color = i < COLOURS.length ? ColorPicker.toString(COLOURS[i]) : "color=black";
            result.put(clusterSizeList.get(i).clusterId, color);
        }

        for (Link link : links) {
            if (link.isSimpleSV()) {
                result.put(link.clusterId(), simpleSvColor(link.type()));
            } else if (link.isLineElement()) {
                result.put(link.clusterId(), toString(LINE));
            }
        }

        return result;
    }

    @NotNull
    private Map<Integer, String> chainMap(@NotNull final List<Link> links) {
        final Map<Integer, String> result = Maps.newHashMap();

        if (!links.isEmpty()) {
            final Link firstLink = links.get(0);
            if (firstLink.isSimpleSV()) {
                final String color = simpleSvColor(firstLink.type());
                links.forEach(x -> result.put(x.chainId(), color));
            } else if(firstLink.isLineElement()) {
                links.forEach(x -> result.put(x.chainId(), toString(LINE)));
            }
            else {
                final List<Integer> chainIds = links.stream().map(Link::chainId).distinct().collect(Collectors.toList());
                for (int i = 0; i < chainIds.size(); i++) {
                    final String color = i < COLOURS.length ? toString(COLOURS[i]) : toString(BLACK);
                    result.put(chainIds.get(i), color);
                }

            }

        }

        return result;
    }

    private static double connectorTransparency(int links) {
        if (links < 10) {
            return 1;
        }

        if (links < 50) {
            return 0.8;
        }

        if (links < 100) {
            return 0.5;
        }

        if (links < 400) {
            return 0.3;
        }

        return 0.1;
    }

}
