package com.hartwig.hmftools.panelbuilder;

import static java.lang.String.format;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CN_BACKBONE_RESOLUTION_MIN;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.Nullable;

// Copy number backbone partition spacing in bases, with optional per-chromosome overrides of the default.
public record CnBackboneResolution(
        int defaultResolution,
        Map<HumanChromosome, Integer> overrides
)
{
    public int forChromosome(final HumanChromosome chromosome)
    {
        return overrides.getOrDefault(chromosome, defaultResolution);
    }

    // defaultKb: global spacing in kb. overridesStr: comma-separated chr:kb pairs, or null for none.
    public static CnBackboneResolution parse(int defaultKb, @Nullable final String overridesStr)
    {
        Map<HumanChromosome, Integer> overrides = (overridesStr == null || overridesStr.isBlank())
                ? Map.of()
                : Arrays.stream(overridesStr.split(","))
                        .map(CnBackboneResolution::parseOverride)
                        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, CnBackboneResolution::rejectDuplicate));
        return new CnBackboneResolution(kbToBases(defaultKb), overrides);
    }

    private static Map.Entry<HumanChromosome, Integer> parseOverride(final String item)
    {
        String[] parts = item.split(":");
        if(parts.length != 2)
        {
            throw new UserInputError("Malformed copy number backbone resolution override: " + item);
        }
        return Map.entry(parseChromosome(parts[0].trim()), kbToBases(parseKb(parts[1].trim())));
    }

    private static HumanChromosome parseChromosome(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
        {
            throw new UserInputError("Invalid chromosome in copy number backbone resolution override: " + chromosome);
        }
        return HumanChromosome.fromString(chromosome);
    }

    private static int parseKb(final String kb)
    {
        try
        {
            return Integer.parseInt(kb);
        }
        catch(NumberFormatException e)
        {
            throw new UserInputError("Invalid resolution in copy number backbone resolution override: " + kb);
        }
    }

    private static int kbToBases(int kb)
    {
        int bases = kb * 1000;
        if(bases < CN_BACKBONE_RESOLUTION_MIN)
        {
            throw new UserInputError(format("Copy number backbone resolution too small: %db (minimum %db)",
                    bases, CN_BACKBONE_RESOLUTION_MIN));
        }
        return bases;
    }

    private static int rejectDuplicate(int a, int b)
    {
        throw new UserInputError("Duplicate chromosome in copy number backbone resolution override");
    }
}
