package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RECALIBRATED_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_SEQ_TECH_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_READ_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_FINAL_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_FINAL_ALT_MAP_QUAL;
import static com.hartwig.hmftools.sage.SageConstants.CORE_AF_FULL_RATIO;
import static com.hartwig.hmftools.sage.SageConstants.CORE_AF_MIN;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LOCAL_PHASE_SET_READ_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.SIMPLE_ALT_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TUMOR_QUALITY_PROB;

import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.evidence.QualCounters;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.ReadSupportCounts;

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

        builder.log10PError((double)round(variant.tumorReadCounters().get(0).logTqp() * 10d) / 10d);

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

        builder.attribute(TUMOR_QUALITY_PROB, primaryRcCounter.tumorQualProbability());
        builder.attribute(MAP_QUAL_FACTOR, primaryRcCounter.mapQualFactor());

        if(variant.nearIndel())
            builder.attribute(NEARBY_INDEL_FLAG, true);

        final VariantContext context = builder.make();
        if(context.isNotFiltered())
        {
            context.getCommonInfo().addFilter(PASS_FILTER);
        }

        return context;
    }

    private static Genotype createGenotype(final ReadContextCounter counter, final String sampleId)
    {
        return createGenotype(counter, sampleId, null);
    }

    public static Genotype createGenotype(final ReadContextCounter counter, final String sampleId, @Nullable final Set<String> existingTags)
    {
        GenotypeBuilder builder = new GenotypeBuilder(sampleId);

        QualCounters qualCounters = counter.qualCounters();

        int depth = counter.depth();
        int altSupport = counter.altSupport();
        int refSupport = counter.refSupport();
        int strongSupport = counter.strongAltSupport();

        int avgMapQuality = depth > 0 ? (int) round(qualCounters.mapQualityTotal() / (double)depth) : 0;
        int avgAltMapQuality = altSupport > 0 ? (int) round(qualCounters.altMapQualityTotal() / (double)altSupport) : 0;
        int avgBaseQuality = depth > 0 ? (int)round(qualCounters.recalibratedBaseQualityTotal() / (double)depth) : 0;
        int avgAltBaseQuality = (int)round(counter.averageAltRecalibratedBaseQuality());

        int avgAltFinalBaseQuality = strongSupport > 0 ? (int)round(qualCounters.altFinalBaseQualityTotal() / (double)strongSupport) : 0;
        int avgAltFinalMapQuality = strongSupport > 0 ? (int)round(qualCounters.altFinalMapQualityTotal() / (double)strongSupport) : 0;

        // NOTE: any field added to the genotype field should also be added to checkGenotypeFields() above, so that if it is only set in
        // append-mode, that it will also get valid default values

        double vaf = counter.vaf();

        // remove core support from AD and AF consideration if it is implausibly high
        ReadSupportCounts readCounts = counter.readCounts();
        double coreProportion = (double)readCounts.Core / readCounts.altSupport();

        if(readCounts.Core >= CORE_AF_MIN && coreProportion >= CORE_AF_FULL_RATIO)
        {
            altSupport -= readCounts.Core;
            vaf = altSupport / (double)readCounts.Total;
        }

        // add back reads with uncertain bases in the core to depth and proportionally to ref and alt support
        int uncertainCoreBaseCount = counter.uncertainCoreBaseCount();

        if(uncertainCoreBaseCount > 0)
        {
            double altRatio = altSupport / (double)depth;
            double refRatio = refSupport / (double)depth;

            depth += uncertainCoreBaseCount;
            refSupport += (int)round(uncertainCoreBaseCount * refRatio);
            altSupport += (int)round(uncertainCoreBaseCount * altRatio);
        }

        builder.DP(depth)
                .AD(new int[] { refSupport, altSupport })
                .attribute(READ_CONTEXT_QUALITY, counter.quality())
                .attribute(READ_CONTEXT_COUNT, readCounts.toArray())
                .attribute(READ_CONTEXT_IMPROPER_PAIR, counter.improperPairCount())
                .attribute(READ_CONTEXT_JITTER, counter.jitter().summary())
                .attribute(AVG_READ_MAP_QUALITY, new int[] { avgMapQuality, avgAltMapQuality })
                .attribute(AVG_RECALIBRATED_BASE_QUAL, new int[] { avgBaseQuality, avgAltBaseQuality })
                .attribute(AVG_SEQ_TECH_BASE_QUAL, (int)counter.averageAltSeqTechBaseQuality())
                .attribute(
                        AVG_EDGE_DISTANCE_PERC, new double[] {
                                counter.readEdgeDistance().avgDistanceFromEdge(),
                                counter.readEdgeDistance().avgAltDistanceFromEdge() })
                .attribute(AVG_FINAL_BASE_QUAL, avgAltFinalBaseQuality)
                .attribute(AVG_FINAL_ALT_MAP_QUAL, avgAltFinalMapQuality)
                .attribute(
                        FRAG_STRAND_BIAS, format("%.3f,%.3f", counter.fragmentStrandBiasNonAlt().bias(), counter.fragmentStrandBiasAlt().bias()))
                .attribute(
                        READ_STRAND_BIAS, format("%.3f,%.3f", counter.readStrandBiasNonAlt().bias(), counter.readStrandBiasAlt().bias()))
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, vaf)
                .attribute(SIMPLE_ALT_COUNT, counter.simpleAltMatches())
                .attribute(MIN_COORDS_COUNT, counter.fragmentCoords().minCount())
                .alleles(NO_CALL);

        if(counter.consensusTypeCounts() != null)
        {
            builder.attribute(UMI_TYPE_COUNTS, counter.consensusTypeCounts());
        }
        else
        {
            builder.attribute(UMI_TYPE_COUNTS, new int[] {counter.depth(), 0, 0, counter.readCounts().strongSupport(), 0, 0});
        }

        if(existingTags != null)
        {
            Set<String> allTags = builder.make().getExtendedAttributes().keySet();

            for(String vcfTag : existingTags)
            {
                if(!allTags.contains(vcfTag))
                    builder.attribute(vcfTag, "");
            }
        }

        return builder.make();
    }
}
