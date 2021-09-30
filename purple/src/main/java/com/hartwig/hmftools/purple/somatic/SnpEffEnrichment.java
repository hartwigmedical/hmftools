package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.GeneNameMapping;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CanonicalAnnotation;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SnpEffEnrichment implements VariantContextEnrichment
{
    private final Consumer<VariantContext> mConsumer;

    private final CanonicalAnnotation mCanonicalAnnotation;
    private final EnsemblDataCache mGeneTransCache;

    public SnpEffEnrichment(
            final Set<String> driverGenes, final EnsemblDataCache geneTransCache, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;

        mGeneTransCache = geneTransCache;
        mGeneTransCache.createGeneNameIdMap();

        Map<String,String> transGeneMap = Maps.newHashMap();

        for(List<GeneData> geneDataList : geneTransCache.getChrGeneDataMap().values())
        {
            for(GeneData geneData : geneDataList)
            {
                List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

                for(TranscriptData tranData : transDataList)
                {
                    transGeneMap.put(tranData.TransName, geneData.GeneName);
                }
            }
        }

        mCanonicalAnnotation = new CanonicalAnnotation(driverGenes, transGeneMap);
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        final VariantImpact variantImpact = formVariantImpact(context);
        VariantImpactSerialiser.writeImpactDetails(context, variantImpact);
        mConsumer.accept(context);
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader header)
    {
        return VariantImpactSerialiser.writeHeader(header);
    }

    private VariantImpact formVariantImpact(@NotNull final VariantContext context)
    {
        boolean phasedInframeIndel = context.isIndel() && context.getAttributeAsInt(SageMetaData.PHASED_INFRAME_INDEL, 0) > 0;

        final List<SnpEffAnnotation> allAnnotations = SnpEffAnnotationParser.fromContext(context);

        final List<SnpEffAnnotation> transcriptAnnotations =
                allAnnotations.stream().filter(SnpEffAnnotation::isTranscriptFeature).collect(Collectors.toList());

        int genesAffected = (int)transcriptAnnotations.stream()
                .map(SnpEffAnnotation::gene)
                .filter(x -> !x.isEmpty())
                .distinct()
                .count();

        String canonicalGeneName = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        boolean canonicalIsSplice = false;
        CodingEffect worstCodingEffect = UNDEFINED;

        if(!transcriptAnnotations.isEmpty())
        {
            final SnpEffAnnotation worstAnnotation = transcriptAnnotations.get(0);
            worstCodingEffect = codingEffect(context, phasedInframeIndel, worstAnnotation);
        }

        final Optional<SnpEffAnnotation> canonicalAnnotation = mCanonicalAnnotation.canonicalSnpEffAnnotation(transcriptAnnotations);
        if(canonicalAnnotation.isPresent())
        {
            final SnpEffAnnotation annotation = canonicalAnnotation.get();
            canonicalGeneName = annotation.gene();
            canonicalEffect = annotation.consequenceString();
            canonicalCodingEffect = codingEffect(context, phasedInframeIndel, annotation);
            canonicalHgvsCodingImpact = annotation.hgvsCoding();
            canonicalHgvsProteinImpact = annotation.hgvsProtein();
            canonicalTranscript = annotation.transcript();
            canonicalIsSplice = canonicalEffect.contains("splice");
        }

        return new VariantImpact(
                canonicalGeneName, canonicalEffect, canonicalTranscript, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, canonicalIsSplice, "", worstCodingEffect, genesAffected);
    }

    private CodingEffect codingEffect(final VariantContext context, boolean phasedInframeIndel, final SnpEffAnnotation annotation)
    {
        // locate the Ensembl data and translate to HmfTranscriptRegion until Pave is complete
        HmfTranscriptRegion transcriptRegion = null;

        GeneData geneData = mGeneTransCache.getGeneDataByName(annotation.gene());

        if(geneData == null)
        {
            // try mapping
            String geneNameNew = mGeneTransCache.getGeneMappings().getNewName(annotation.gene());

            if(geneNameNew != null)
                geneData = mGeneTransCache.getGeneDataByName(geneNameNew);
        }

        if(geneData != null)
        {
            List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

            if(transDataList == null)
            {
                TranscriptData transData = transDataList.stream()
                        .filter(x -> x.TransName.equals(annotation.featureID())).findFirst().orElse(null);

                if(transData != null)
                    transcriptRegion = HmfTranscriptRegionUtils.fromTranscript(geneData, transData);
            }
        }

        CodingEffect effect = CodingEffectFactory.effect(context, transcriptRegion, annotation.consequences());
        return phasedInframeIndel && effect.equals(CodingEffect.NONSENSE_OR_FRAMESHIFT) ? CodingEffect.MISSENSE : effect;
    }

    @Override
    public void flush() {}

}
