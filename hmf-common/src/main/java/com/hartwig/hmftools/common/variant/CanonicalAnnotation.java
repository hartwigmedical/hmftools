package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;

import org.jetbrains.annotations.NotNull;

public class CanonicalAnnotation
{
    private final Set<String> mDriverCatalogGenes;
    private final Map<String, String> mTranscriptGeneMap;

    /*
    public CanonicalAnnotation(final Set<String> driverGenes, final List<HmfTranscriptRegion> transcripts)
    {
        mDriverCatalogGenes = driverGenes;

        mTranscriptGeneMap = transcripts.stream()
                .collect(Collectors.toMap(TranscriptRegion::transName, TranscriptRegion::geneName));
    }
    */

    public CanonicalAnnotation(final Set<String> driverGenes, final Map<String,String> transGeneMap)
    {
        mDriverCatalogGenes = driverGenes;
        mTranscriptGeneMap = transGeneMap;
    }

    public Optional<SnpEffAnnotation> canonicalSnpEffAnnotation(final List<SnpEffAnnotation> allAnnotations)
    {
        final List<SnpEffAnnotation> transcriptAnnotations = allAnnotations.stream()
                .filter(SnpEffAnnotation::isTranscriptFeature)
                .collect(Collectors.toList());

        return pickCanonicalFavourDriverGene(transcriptAnnotations);
    }

    @VisibleForTesting
    @NotNull
    <T extends SnpEffAnnotation> Optional<T> pickCanonicalFavourDriverGene(List<T> annotations)
    {
        List<T> canonicalAnnotations = annotations.stream()
                .filter(annotation -> mTranscriptGeneMap.containsKey(trimEnsembleVersion(annotation.transcript())))
                .collect(Collectors.toList());

        if(!canonicalAnnotations.isEmpty())
        {
            Optional<T> canonicalOnDriverGene =
                    canonicalAnnotations.stream().filter(annotation -> mDriverCatalogGenes.contains(annotation.gene())).findFirst();
            if(canonicalOnDriverGene.isPresent())
            {
                return canonicalOnDriverGene;
            }

            return Optional.of(canonicalAnnotations.get(0));
        }

        return Optional.empty();
    }

    static String trimEnsembleVersion(final String transcriptId)
    {
        if(transcriptId.startsWith("EN") && transcriptId.contains("."))
        {
            return transcriptId.substring(0, transcriptId.indexOf("."));
        }

        return transcriptId;
    }
}
