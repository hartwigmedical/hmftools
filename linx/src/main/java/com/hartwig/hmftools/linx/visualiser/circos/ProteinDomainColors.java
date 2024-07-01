package com.hartwig.hmftools.linx.visualiser.circos;

import java.awt.Color;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;

import org.jetbrains.annotations.NotNull;

public class ProteinDomainColors
{
    private static final float SAT = 0.5f;
    private static final float BRIGHT = 0.8392f;
    private static final Color UTR = new Color(128, 128, 128);
    private static final Color UTR_DOWN = new Color(107, 174, 214);
    private static final Color UTR_UP = new Color(214, 144, 107);

    private final Map<String, Color> proteinColorMap = Maps.newLinkedHashMap();

    public ProteinDomainColors(final List<VisProteinDomain> proteinDomains)
    {
        this(proteinDomains.stream().map(VisProteinDomain::name).collect(Collectors.toSet()));
    }

    ProteinDomainColors(final Set<String> proteinDomains)
    {
        proteinColorMap.put(VisProteinDomain.PROTEIN_DOMAIN_UTR, UTR);

        final List<String> newProteinDomains =
                proteinDomains.stream().filter(x -> !proteinColorMap.containsKey(x)).collect(Collectors.toList());

        for(int i = 0; i < newProteinDomains.size(); i++)
        {
            proteinColorMap.put(newProteinDomains.get(i), getFloatingColor(i, newProteinDomains.size()));
        }

    }

    @NotNull
    public Color color(final VisProteinDomain proteinDomain)
    {
        if(proteinDomain.name().equals(VisProteinDomain.PROTEIN_DOMAIN_UTR))
        {
            return proteinDomain.start() == 1 ? UTR_DOWN : UTR_UP;
        }

        return proteinColorMap.get(proteinDomain.name());
    }

    private static Color getFloatingColor(int i, int maxDomains)
    {
        return Color.getHSBColor(hue(i, maxDomains), SAT, BRIGHT);
    }

    static float hue(int i, int maxDomains)
    {
        float colorDistance = 1f / (float) Math.ceil((maxDomains + 2) / 2f) / 2f;
        return i % 2 == 0
                ? 21f / 360f + (i / 2 + 1) * colorDistance
                : 201f / 360f + (float) Math.ceil(i / 2f) * colorDistance;
    }

    @Override
    @NotNull
    public String toString()
    {
        final StringBuilder builder = new StringBuilder();
        final List<Color> colors = Lists.newArrayList(proteinColorMap.values());
        for(int i = 0; i < colors.size(); i++)
        {
            final Color color = colors.get(i);
            builder.append("private static final Color COLOR").append(i).append(" = new Color(")
                    .append(color.getRed()).append(",")
                    .append(color.getGreen()).append(",")
                    .append(color.getBlue())
                    .append(");").append("\n");
        }

        return builder.toString();
    }

}
