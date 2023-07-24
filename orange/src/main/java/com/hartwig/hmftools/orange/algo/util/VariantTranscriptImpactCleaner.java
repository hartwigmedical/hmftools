package com.hartwig.hmftools.orange.algo.util;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VariantTranscriptImpactCleaner
{
    private VariantTranscriptImpactCleaner()
    {
    }

    /**
     * When VariantTranscriptImpacts are created from the VCF file, it sometimes (incorrectly) parses the square array brackets.
     * These parsed brackets are then included in the impacts fields and here we remove them.
     * TODO maybe fix this bug upstream in the VariantTranscriptImpact class?
     */
    @NotNull
    public static VariantTranscriptImpact cleanFields(@NotNull VariantTranscriptImpact impact)
    {
        var cleanedGeneId = stripSquareBracketsAndWhiteSpace(impact.GeneId);
        var cleanedGeneName = stripSquareBracketsAndWhiteSpace(impact.GeneName);
        var cleanedTranscript = stripSquareBracketsAndWhiteSpace(impact.Transcript);
        var cleanedHgvsCoding = stripSquareBracketsAndWhiteSpace(impact.HgvsCoding);
        var cleanedHgvsProtein = stripSquareBracketsAndWhiteSpace(impact.HgvsProtein);
        var cleanedEffects = stripSquareBracketsAndWhiteSpace(impact.Effects);
        return new VariantTranscriptImpact(
                cleanedGeneId,
                cleanedGeneName,
                cleanedTranscript,
                cleanedEffects,
                impact.SpliceRegion,
                cleanedHgvsCoding,
                cleanedHgvsProtein);
    }

    @Nullable
    private static String stripSquareBracketsAndWhiteSpace(String string)
    {
        if(string == null)
        {
            return null;
        }
        string = string.startsWith("[") ? string.substring(1) : string;
        string = string.endsWith("]") ? string.substring(0, string.length() - 1) : string;
        return string.strip();
    }
}
