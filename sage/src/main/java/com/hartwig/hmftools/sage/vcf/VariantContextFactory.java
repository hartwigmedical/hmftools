package com.hartwig.hmftools.sage.vcf;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_ALT_MAP_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LOCAL_PHASE_SET_READ_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.QUAL_MODEL_TYPE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.SIMPLE_ALT_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TUMOR_QUALITY_PROB;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.evidence.QualCounters;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.common.SageVariant;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public final class VariantContextFactory
{
    public static final List<Allele> NO_CALL = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    public static VariantContext create(final SageVariant variant, final List<String> referenceIds, final List<String> tumorIds)
    {
        try
        {
            final List<Genotype> genotypes = Lists.newArrayList();
            for(int i = 0; i < variant.referenceReadCounters().size(); i++)
            {
                ReadContextCounter refCounter = variant.referenceReadCounters().get(i);
                genotypes.add(createGenotype(refCounter, referenceIds.get(i)));
            }

            for(int i = 0; i < variant.tumorReadCounters().size(); i++)
            {
                ReadContextCounter tumorCounter = variant.tumorReadCounters().get(i);
                genotypes.add(createGenotype(tumorCounter, tumorIds.get(i)));
            }

            return createContext(variant, genotypes);
        }
        catch(Exception e)
        {
            SG_LOGGER.error("var({}) failed to create VCF context: {}", variant.candidate().readContext(), e.toString());
            System.exit(1);
            return null;
        }
    }

    private static VariantContext createContext(final SageVariant variant, final List<Genotype> genotypes)
    {
        VariantContextBuilder builder = CandidateSerialisation.toContext(variant.candidate());

        builder.log10PError(variant.totalQuality() / -10d);
        builder.genotypes(genotypes);
        builder.filters(variant.filtersStringSet());

        if(variant.hasLocalPhaseSets())
        {
            // write in order - PAVE will use the first only, but all are loaded to the DB
            builder.attribute(LOCAL_PHASE_SET, variant.localPhaseSets());

            List<Integer> readCountsRaw = variant.localPhaseSetCounts();
            List<Integer> readCountTotals = Lists.newArrayListWithExpectedSize(readCountsRaw.size());
            readCountsRaw.forEach(x -> readCountTotals.add(x));
            builder.attribute(LOCAL_PHASE_SET_READ_COUNT, readCountTotals);
        }

        if(variant.mixedGermlineImpact() > 0)
        {
            builder.attribute(MIXED_SOMATIC_GERMLINE, variant.mixedGermlineImpact());
        }

        ReadContextCounter primaryRcCounter = variant.tumorReadCounters().get(0);

        builder.attribute(MAX_READ_EDGE_DISTANCE, primaryRcCounter.readEdgeDistance().maxAltDistanceFromEdge());

        builder.attribute(
                AVG_READ_EDGE_DISTANCE,
                new int[] { primaryRcCounter.readEdgeDistance().avgDistanceFromEdge(),
                        primaryRcCounter.readEdgeDistance().avgAltDistanceFromEdge() } );

        builder.attribute(TUMOR_QUALITY_PROB, primaryRcCounter.tumorQualProbability());
        builder.attribute(MAP_QUAL_FACTOR, primaryRcCounter.mapQualFactor());

        if(primaryRcCounter.ultimaQualModel() != null)
        {
            builder.attribute(QUAL_MODEL_TYPE, primaryRcCounter.ultimaQualModel().type().toString());
        }

        final VariantContext context = builder.make();
        if(context.isNotFiltered())
        {
            context.getCommonInfo().addFilter(PASS);
        }

        return context;
    }

    public static Genotype createGenotype(final ReadContextCounter counter, final String sampleId)
    {
        GenotypeBuilder builder = new GenotypeBuilder(sampleId);

        QualCounters qualCounters = counter.qualCounters();

        int depth = counter.depth();
        int altSupport = counter.altSupport();
        int strongSupport = counter.strongAltSupport();

        int avgMapQuality = depth > 0 ? (int) round(qualCounters.mapQualityTotal() / (double)depth) : 0;
        int avgAltMapQuality = altSupport > 0 ? (int) round(qualCounters.altMapQualityTotal() / (double)altSupport) : 0;
        int avgBaseQuality = depth > 0 ? (int)round(qualCounters.baseQualityTotal() / (double)depth) : 0;
        int avgAltBaseQuality = (int)round(counter.averageAltBaseQuality());

        int avgAltModifiedBaseQuality = strongSupport > 0 ? (int)round(qualCounters.modifiedAltBaseQualityTotal() / (double)strongSupport) : 0;
        int avgAltModifiedMapQuality = strongSupport > 0 ? (int)round(qualCounters.altModifiedMapQualityTotal() / (double)strongSupport) : 0;

        builder.DP(depth)
                .AD(new int[] { counter.refSupport(), altSupport })
                .attribute(READ_CONTEXT_QUALITY, counter.quality())
                .attribute(READ_CONTEXT_COUNT, counter.counts())
                .attribute(READ_CONTEXT_IMPROPER_PAIR, counter.improperPairCount())
                .attribute(READ_CONTEXT_JITTER, counter.jitter().summary())
                .attribute(AVG_MAP_QUALITY, new int[] { avgMapQuality, avgAltMapQuality })
                .attribute(AVG_BASE_QUAL, new int[] { avgBaseQuality, avgAltBaseQuality })
                .attribute(AVG_MODIFIED_BASE_QUAL, avgAltModifiedBaseQuality)
                .attribute(AVG_MODIFIED_ALT_MAP_QUAL, avgAltModifiedMapQuality)
                .attribute(
                        FRAG_STRAND_BIAS, format("%.3f,%.3f", counter.fragmentStrandBiasNonAlt().bias(), counter.fragmentStrandBiasAlt().bias()))
                .attribute(
                        READ_STRAND_BIAS, format("%.3f,%.3f", counter.readStrandBiasNonAlt().bias(), counter.readStrandBiasAlt().bias()))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, counter.vaf())
                .attribute(SIMPLE_ALT_COUNT, counter.simpleAltMatches())
                .alleles(NO_CALL);

        if(counter.umiTypeCounts() != null)
        {
            builder.attribute(UMI_TYPE_COUNTS, counter.umiTypeCounts());
        }

        return builder.make();
    }
}
