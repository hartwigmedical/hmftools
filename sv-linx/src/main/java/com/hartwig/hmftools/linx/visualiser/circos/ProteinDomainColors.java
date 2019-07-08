package com.hartwig.hmftools.linx.visualiser.circos;

import java.awt.Color;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomains;

import org.jetbrains.annotations.NotNull;

public class ProteinDomainColors
{

    // Top 10 were selected with the following query:
    //
    // proteinDomain = read.table(file = "~/hmf/analysis/fusions/SVA_VIS_PROTEIN_DOMAINS.tsv", sep = "\t", header = T) %>%
    //     select(Transcript, Info) %>% distinct() %>%
    //     group_by(Info) %>% count() %>% ungroup() %>%
    //     top_n(10, n) %>% arrange(-n)

    private static final float FIXED_HUE = 0.166f;
    private static final float SAT = 0.7f;
    private static final float BRIGHT = 0.8f;
    private static final Color UTR = new Color(128,128,128);

    private final Map<String, Color> proteinColorMap = Maps.newLinkedHashMap();

    ProteinDomainColors(@NotNull final Set<String> proteinDomains)
    {
        proteinColorMap.put("Protein kinase domain", getFixedColor(0));
        proteinColorMap.put("Immunoglobulin-like domain", getFixedColor(1));
        proteinColorMap.put("Zinc finger; C2H2", getFixedColor(2));
        proteinColorMap.put("Zinc finger; RING-type", getFixedColor(3));
        proteinColorMap.put("Pleckstrin homology domain", getFixedColor(4));
        proteinColorMap.put("RNA recognition motif domain", getFixedColor(5));
        proteinColorMap.put("SH2 domain", getFixedColor(6));
        proteinColorMap.put("Ankyrin repeat", getFixedColor(7));
        proteinColorMap.put("Ankyrin repeat-containing domain", getFixedColor(8));
        proteinColorMap.put("Fibronectin type III", getFixedColor(9));
        proteinColorMap.put(ProteinDomains.UTR, UTR);

        final List<String> newProteinDomains =
                proteinDomains.stream().filter(x -> !proteinColorMap.keySet().contains(x)).collect(Collectors.toList());

        for (int i = 0; i < newProteinDomains.size(); i++)
        {
            proteinColorMap.put(newProteinDomains.get(i), getFloatingColor(i, newProteinDomains.size()));
        }

    }

    @NotNull
    public String rgb(@NotNull final String proteinDomain)
    {
        final Color color = color(proteinDomain);
        return "(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ",0.8)";
    }

    @NotNull
    public Color color(@NotNull final String proteinDomain)
    {
        return proteinColorMap.get(proteinDomain);
    }

    private static Color getFloatingColor(int i, int maxDomains)
    {
        return Color.getHSBColor(FIXED_HUE + ((i + 1) * (1f - FIXED_HUE)) / (maxDomains + 1), SAT, BRIGHT);
    }

    private static Color getFixedColor(int i)
    {
        return Color.getHSBColor(i * FIXED_HUE / 9f, SAT, BRIGHT);
    }

    @Override
    @NotNull
    public String toString()
    {
        final StringBuilder builder = new StringBuilder();
        final List<Color> colors = Lists.newArrayList(proteinColorMap.values());
        for (int i = 0; i < colors.size(); i++)
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
