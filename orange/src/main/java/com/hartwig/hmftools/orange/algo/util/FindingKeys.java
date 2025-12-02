package com.hartwig.hmftools.orange.algo.util;

import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FindingKeys
{
    public enum SampleType
    {
        GERMLINE, SOMATIC
    }

    public static String smallVariant(@NotNull SampleType sampleType, @NotNull PurpleVariant variant,
            @NotNull PurpleTranscriptImpact transcriptImpact, boolean isCanonical)
    {
        return String.format("smallVariant[%s %s %s]", sampleType, geneTranscriptLabel(variant.gene(), isCanonical,
                transcriptImpact.transcript()), impact(transcriptImpact));
    }

    public static String gainDeletion(@NotNull SampleType sampleType, String gene, CopyNumberInterpretation copyNumberInterpretation,
            boolean isCanonical, String transcriptId)
    {
        return String.format("gainDeletion[%s %s %s]", sampleType, geneTranscriptLabel(gene, isCanonical, transcriptId),
                copyNumberInterpretation.name());
    }

    public static String disruption(@NotNull SampleType sampleType, LinxBreakend breakend)
    {
        return String.format("disruption[%s %s %d]",
                sampleType,
                geneTranscriptLabel(breakend.gene(), breakend.isCanonical(), breakend.transcript()),
                breakend.svId());
    }

    public static String fusion(@NotNull SampleType sampleType, LinxFusion fusion)
    {
        return String.format("fusion[%s %s %s]", sampleType, fusion.geneStart(), fusion.geneEnd());
    }

    public static String virus(VirusInterpreterEntry virus) {
        String label = virus.interpretation() + (virus.interpretation() == VirusInterpretation.HPV ? " (" + virus.name() + ")" : "");
        return String.format("virus[%s]", label);
    }

    public static String microsatelliteStability(MicrosatelliteStatus status) {
        return String.format("microsatelliteStability[%s]", status.name());
    }

    public static String homologousRecombination(ChordStatus status) {
        return String.format("homologousRecombination[%s]", status.name());
    }

    public static String tumorMutationStatus(PurpleTumorMutationalStatus tmbStatus, PurpleTumorMutationalStatus tmlStatus) {
        return String.format("tumorMutationStatus[TMB_%s TML_%s]", tmbStatus, tmlStatus);
    }

    public static String predictedTumorOrigin(@NotNull String cancerType) {
        return String.format("predictedTumorOrigin[%s]", cancerType);
    }

    // only show transcript ID for non canonical transcripts
    private static String geneTranscriptLabel(String gene, boolean isCanonical, String transcriptId) {
        return isCanonical ? gene : String.format("%s(%s)", gene, transcriptId);
    }

    private static String impact(@NotNull PurpleTranscriptImpact transcriptImpact) {
        return determineVariantAnnotation(transcriptImpact.hgvsCodingImpact(),
                transcriptImpact.hgvsProteinImpact(),
                transcriptImpact.effects().stream().map(Enum::toString).collect(Collectors.joining("&")),
                transcriptImpact.codingEffect() == PurpleCodingEffect.SPLICE,
                transcriptImpact.effects().contains(PurpleVariantEffect.UPSTREAM_GENE));
    }

    @NotNull
    public static String determineVariantAnnotation(@Nullable String hgvsCoding, @Nullable String hgvsProtein, @NotNull String effects,
            boolean isSplice, boolean isUpstream) {
        if (!Strings.isNullOrEmpty(hgvsProtein) && !hgvsProtein.equals("p.?")) {
            return hgvsProtein;
        }

        if (!Strings.isNullOrEmpty(hgvsCoding)) {
            return isSplice ? hgvsCoding + " splice" : hgvsCoding;
        }

        return isUpstream ? "upstream" : effects;
    }
}
