package com.hartwig.hmftools.sage.vcf;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Doubles.round;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.AVG_NM_COUNT;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.LOCAL_PHASE_SET_READ_COUNT;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.OLD_INDEL_DEDUP_FLAG;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_ALLELIC_BASE_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_ALLELIC_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.RAW_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.READ_STRAND_BIAS;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.sage.append.CandidateSerialization;
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
        final List<Genotype> genotypes = Lists.newArrayList();
        for(int i = 0; i < variant.normalReadCounters().size(); i++)
        {
            ReadContextCounter normalContext = variant.normalReadCounters().get(i);
            genotypes.add(createGenotype(normalContext, referenceIds.get(i)));
        }

        for(int i = 0; i < variant.tumorReadCounters().size(); i++)
        {
            ReadContextCounter tumorContext = variant.tumorReadCounters().get(i);
            genotypes.add(createGenotype(tumorContext, tumorIds.get(i)));
        }

        return createContext(variant, genotypes);
    }

    private static VariantContext createContext(final SageVariant variant, final List<Genotype> genotypes)
    {
        final VariantContextBuilder builder = CandidateSerialization.toContext(variant.candidate())
                .log10PError(variant.totalQuality() / -10d)
                .genotypes(genotypes)
                .filters(variant.filters());

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

        if(variant.dedupIndelDiff())
        {
            builder.attribute(OLD_INDEL_DEDUP_FLAG, true);
        }

        builder.attribute(MAX_READ_EDGE_DISTANCE, variant.tumorReadCounters().get(0).readEdgeDistance().maxAltDistanceFromUnclippedEdge());

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

        int depth = counter.depth();
        int altSupport = counter.altSupport();

        int avgMapQuality = depth > 0 ? (int)Math.round(counter.mapQualityTotal() / (double)depth) : 0;
        int avgAltMapQuality = altSupport > 0 ? (int)Math.round(counter.altMapQualityTotal() / (double)altSupport) : 0;

        double avgNmCount = depth > 0 ? round(counter.nmCountTotal() / (double)depth, 1) : 0;
        double avgAltNmCount = altSupport > 0 ? round(counter.altNmCountTotal() / (double)altSupport, 1) : 0;

        builder.DP(depth)
                .AD(new int[] { counter.refSupport(), altSupport })
                .attribute(READ_CONTEXT_QUALITY, counter.quality())
                .attribute(READ_CONTEXT_COUNT, counter.counts())
                .attribute(READ_CONTEXT_IMPROPER_PAIR, counter.improperPairCount())
                .attribute(READ_CONTEXT_JITTER, counter.jitter())
                .attribute(AVG_MAP_QUALITY, new int[] { avgMapQuality, avgAltMapQuality })
                .attribute(AVG_NM_COUNT, new double[] { avgNmCount, avgAltNmCount })
                .attribute(RAW_ALLELIC_DEPTH, new int[] { counter.rawRefSupport(), counter.rawAltSupport() })
                .attribute(RAW_ALLELIC_BASE_QUALITY, new int[] { counter.rawRefBaseQualityTotal(), counter.rawAltBaseQualityTotal() })
                .attribute(RAW_DEPTH, counter.rawDepth())
                .attribute(
                        FRAG_STRAND_BIAS, format("%.3f,%.3f", counter.fragmentStrandBiasAlt().bias(), counter.fragmentStrandBiasRef().bias()))
                .attribute(
                        READ_STRAND_BIAS, format("%.3f,%.3f", counter.readStrandBiasAlt().bias(), counter.readStrandBiasRef().bias()))
                .attribute(AVG_BASE_QUAL, (int)counter.averageAltBaseQuality())
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, counter.vaf())
                .alleles(NO_CALL);

        if(counter.umiTypeCounts() != null)
        {
            builder.attribute(UMI_TYPE_COUNTS, counter.umiTypeCounts());
        }

        return builder.make();
    }
}
