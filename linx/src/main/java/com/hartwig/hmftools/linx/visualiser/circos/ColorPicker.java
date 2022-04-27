package com.hartwig.hmftools.linx.visualiser.circos;

import static java.util.stream.Collectors.toMap;

import java.awt.Color;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.jetbrains.annotations.NotNull;

public class ColorPicker
{
    private static final Map<String, Color> COLOR_MAP = Maps.newHashMap();
    static {
        COLOR_MAP.put("1", new Color(128, 125, 186));
        COLOR_MAP.put("2", new Color(145, 142, 179));
        COLOR_MAP.put("3", new Color(161, 159, 173));
        COLOR_MAP.put("4", new Color(179, 176, 166));
        COLOR_MAP.put("5", new Color(196, 193, 160));
        COLOR_MAP.put("6", new Color(213, 210, 153));

        COLOR_MAP.put("7", new Color(230, 228, 147));
        COLOR_MAP.put("8", new Color(202, 218, 138));
        COLOR_MAP.put("9", new Color(175, 209, 129));
        COLOR_MAP.put("10", new Color(147, 199, 120));
        COLOR_MAP.put("11", new Color(120, 190, 111));
        COLOR_MAP.put("12", new Color(92, 180, 102));

        COLOR_MAP.put("13", new Color(65, 171, 93));
        COLOR_MAP.put("14", new Color(65, 166, 110));
        COLOR_MAP.put("15", new Color(65, 162, 128));
        COLOR_MAP.put("16", new Color(65, 158, 145));
        COLOR_MAP.put("17", new Color(65, 154, 163));
        COLOR_MAP.put("18", new Color(65, 150, 180));

        COLOR_MAP.put("19", new Color(66, 146, 198));
        COLOR_MAP.put("20", new Color(76, 142, 196));
        COLOR_MAP.put("21", new Color(86, 139, 194));
        COLOR_MAP.put("22", new Color(97, 135, 192));
        COLOR_MAP.put("X", new Color(107, 132, 190));
        COLOR_MAP.put("Y", new Color(117, 128, 188));
    }

    public static Color contigColour(@NotNull final String contig)
    {
        return COLOR_MAP.getOrDefault(contig.replace("chr", ""), Color.BLACK);
    }

    @NotNull
    public static String hexContigColor(@NotNull final String contig)
    {
        return hexColor(contigColour(contig));
    }

    @NotNull
    public static String hexColor(@NotNull final Color color)
    {
        return String.format("#%02X%02X%02X", color.getRed(), color.getGreen(), color.getBlue());
    }


    private static final Color BLACK = new Color(5, 5, 5);

    private static final Color COLOR1 = new Color(106, 61, 154);
    private static final Color COLOR2 = new Color(140, 81, 10);
    private static final Color COLOR3 = new Color(1, 102, 94);
    private static final Color COLOR4 = new Color(255, 127, 0);
    private static final Color COLOR5 = new Color(212, 193, 23);
    private static final Color COLOR6 = new Color(31, 120, 180);
    private static final Color COLOR7 = new Color(51, 160, 44);
    private static final Color COLOR8 = new Color(152, 51, 160);

    private static final Color DEL = new Color(251, 154, 153);
    private static final Color DUP = new Color(178, 223, 138);
    private static final Color INS = new Color(255, 255, 153);
    private static final Color LINE = new Color(97, 171, 227);
    private static final Color DOUBLE_MINUTE = new Color(227, 26, 28);

    private static final Color[] COLOURS = new Color[] { COLOR1, COLOR2, COLOR3, COLOR4, COLOR5, COLOR6, COLOR7, COLOR8 };

    private final boolean clusterMode;
    private final Map<Integer, String> colorMap;
    private final double connectorTransparency;

    @NotNull
    public static ColorPicker clusterColors(@NotNull final List<VisSvData> links)
    {
        return new ColorPicker(colorsByCluster(links), true, connectorTransparency(links.size()));
    }

    @NotNull
    public static ColorPicker chainColors(@NotNull final List<VisSvData> links)
    {
        return new ColorPicker(colorsByChain(links), false, connectorTransparency(links.size()));
    }

    private ColorPicker(@NotNull final Map<Integer, String> colorMap, final boolean clusterMode, final double connectorTransparency)
    {
        this.clusterMode = clusterMode;
        this.colorMap = colorMap;
        this.connectorTransparency = connectorTransparency;
    }

    @NotNull
    public String transparentColor(final int clusterId, final int chainId)
    {
        String opaqueColor = color(clusterId, chainId);
        return opaqueColor.replace(")", "," + connectorTransparency + ")");
    }

    @NotNull
    public String color(final int clusterId, final int chainId)
    {
        if (clusterId == -1 || chainId == -1)
        {
            return toString(BLACK);
        }

        return clusterMode ? colorMap.get(clusterId) : colorMap.get(chainId);
    }

    @NotNull
    private static String toString(@NotNull final Color color)
    {
        return "color=(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ")";
    }

    @NotNull
    private static String simpleSvColor(@NotNull final StructuralVariantType type)
    {
        switch (type)
        {
            case DEL:
                return toString(DEL);
            case DUP:
                return toString(DUP);
        }

        return toString(INS);
    }

    private static class ClusterSize
    {

        final int clusterId;
        final long count;

        ClusterSize(final int clusterId, final long count)
        {
            this.clusterId = clusterId;
            this.count = count;
        }
    }

    @NotNull
    private static Map<Integer, String> colorsByCluster(@NotNull final List<VisSvData> links)
    {
        final Map<Integer, String> result = Maps.newHashMap();

        final Comparator<ClusterSize> longComparator = Comparator.<ClusterSize>comparingLong(x -> x.count).reversed();

        final List<ClusterSize> clusterSizeList = Lists.newArrayList();
        links.stream()
                .map(x -> x.ClusterId)
                .collect(toMap(x -> x, (x) -> 1L, Math::addExact))
                .forEach((key, value) -> clusterSizeList.add(new ClusterSize(key, value)));

        clusterSizeList.sort(longComparator);

        for (int i = 0; i < clusterSizeList.size(); i++)
        {
            String color = i < COLOURS.length ? ColorPicker.toString(COLOURS[i]) : "color=black";
            result.put(clusterSizeList.get(i).clusterId, color);
        }

        for (VisSvData link : links)
        {
            if (link.isSimpleSV())
            {
                result.put(link.ClusterId, simpleSvColor(link.Type));
            }
            else if (link.isLineElement())
            {
                result.put(link.ClusterId, toString(LINE));
            }
        }

        return result;
    }

    @NotNull
    private static Map<Integer, String> colorsByChain(@NotNull final List<VisSvData> links)
    {
        final Map<Integer,String> result = Maps.newHashMap();

        if (!links.isEmpty())
        {
            final VisSvData firstLink = links.get(0);

            if(firstLink.isSimpleSV() && !firstLink.InDoubleMinute)
            {
                final String color = simpleSvColor(firstLink.Type);
                links.forEach(x -> result.put(x.ChainId, color));
            }
            else if(firstLink.isLineElement())
            {
                links.forEach(x -> result.put(x.ChainId, toString(LINE)));
            }
            else
            {
                final List<Integer> dmChainIds = links.stream().filter(x -> x.InDoubleMinute).map(x -> x.ChainId)
                        .distinct().collect(Collectors.toList());

                final List<Integer> chainIds = links.stream().map(x -> x.ChainId).distinct().collect(Collectors.toList());

                for(int i = 0; i < chainIds.size(); i++)
                {
                    int chainId = chainIds.get(i);
                    if(dmChainIds.contains(chainId))
                    {
                        result.put(chainId, toString(DOUBLE_MINUTE));
                    }
                    else
                    {
                        final String color = i < COLOURS.length ? toString(COLOURS[i]) : toString(BLACK);
                        result.put(chainId, color);
                    }
                }
            }
        }

        return result;
    }

    private static double connectorTransparency(int links)
    {
        if (links < 10)
        {
            return 1;
        }

        if (links < 50)
        {
            return 0.8;
        }

        if (links < 100)
        {
            return 0.6;
        }

        if (links < 200)
        {
            return 0.4;
        }

        if (links < 400)
        {
            return 0.2;
        }

        return 0.1;
    }

}
