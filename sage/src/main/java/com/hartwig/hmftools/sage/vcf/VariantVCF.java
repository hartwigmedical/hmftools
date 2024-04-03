package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_INDEL_FILTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_MATCH;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_MIXED_GERMLINE_SOMATIC_FILTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_MNV_FILTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_SNV_MNV_FILTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LOCAL_PHASE_SET_READ_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LPS_READ_COUNT_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.QUAL_MODEL_TYPE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.QUAL_MODEL_TYPE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.RAW_SUPPORT_BASE_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.RAW_SUPPORT_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VcfTags.RAW_DEPTH;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_AF_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INDEX;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_MICRO_HOMOLOGY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_MICRO_HOMOLOGY_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_REPEAT_SEQUENCE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_RIGHT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TOTAL_RAW_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TOTAL_RAW_BASE_QUAL_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.VERSION_META_DATA;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.variant.SageVcfTags;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.filter.SoftFilter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VariantVCF implements AutoCloseable
{
    private final VariantContextWriter mWriter;
    private final SomaticRefContextEnrichment mRefContextEnrichment;

    public VariantVCF(
            final IndexedFastaSequenceFile reference, final SageConfig config,
            final List<String> tumorIds, final List<String> referenceIds)
    {
        final SAMSequenceDictionary sequenceDictionary = reference.getSequenceDictionary();

        mWriter = new VariantContextWriterBuilder().setOutputFile(config.OutputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        mRefContextEnrichment = new SomaticRefContextEnrichment(reference);

        final List<String> samples = Lists.newArrayList();
        samples.addAll(referenceIds);
        samples.addAll(tumorIds);

        VCFHeader vcfHeader = createHeader(config.Version, samples, config.Sequencing.Type == SequencingType.ULTIMA);

        final SAMSequenceDictionary condensedDictionary = new SAMSequenceDictionary();
        for(SAMSequenceRecord sequence : sequenceDictionary.getSequences())
        {
            if(HumanChromosome.contains(sequence.getContig()) || MitochondrialChromosome.contains(sequence.getContig()))
            {
                condensedDictionary.addSequence(sequence);
            }
        }

        vcfHeader.setSequenceDictionary(condensedDictionary);
        mWriter.writeHeader(vcfHeader);
    }

    public VariantVCF(final IndexedFastaSequenceFile reference, final SageConfig config, final VCFHeader existingHeader)
    {
        Set<VCFHeaderLine> headerLines = existingHeader.getMetaDataInInputOrder();
        List<String> samples = Lists.newArrayList(existingHeader.getGenotypeSamples());
        samples.addAll(config.ReferenceIds);

        final VCFHeader newHeader = new VCFHeader(headerLines, samples);

        mWriter = new VariantContextWriterBuilder().setOutputFile(config.OutputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(reference.getSequenceDictionary())
                .build();
        mWriter.writeHeader(newHeader);

        mRefContextEnrichment = new SomaticRefContextEnrichment(reference);
    }

    public void write(final VariantContext context)
    {
        mRefContextEnrichment.processVariant(context);
        mWriter.add(context);
    }

    public static VCFHeader createHeader(final String version, final List<String> allSamples, final boolean includeModelData)
    {
        VCFHeader header = new VCFHeader(Collections.emptySet(), allSamples);

        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, LOCAL_PHASE_SET_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_COUNT, 1, VCFHeaderLineType.Integer, READ_CONTEXT_REPEAT_COUNT_DESC));

        // standard fields
        header.addMetaDataLine(new VCFHeaderLine(VERSION_META_DATA, version));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));

        // info fields
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_EVENTS, 1, VCFHeaderLineType.Integer, READ_CONTEXT_EVENTS_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, "Read context core"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_INDEX, 1, VCFHeaderLineType.Integer, "Read context index"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_LEFT_FLANK, 1, VCFHeaderLineType.String, "Read context left flank"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_RIGHT_FLANK, 1, VCFHeaderLineType.String, "Read context right flank"));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_SEQUENCE, 1, VCFHeaderLineType.String, READ_CONTEXT_REPEAT_SEQUENCE_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_MICRO_HOMOLOGY, 1, VCFHeaderLineType.String, READ_CONTEXT_MICRO_HOMOLOGY_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                MIXED_SOMATIC_GERMLINE, 1, VCFHeaderLineType.Integer, MIXED_SOMATIC_GERMLINE_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET_READ_COUNT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, LPS_READ_COUNT_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(MAX_READ_EDGE_DISTANCE, 1, VCFHeaderLineType.Integer, MAX_READ_EDGE_DISTANCE_DESC));

        if(includeModelData)
        {
            header.addMetaDataLine(new VCFInfoHeaderLine(QUAL_MODEL_TYPE, 1, VCFHeaderLineType.String, QUAL_MODEL_TYPE_DESC));
        }

        SageVcfTags.addRefContextHeader(header);

        // genotype fields
        header.addMetaDataLine(new VCFFormatHeaderLine(
                VCFConstants.ALLELE_FREQUENCY_KEY, 1, VCFHeaderLineType.Float, READ_CONTEXT_AF_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_JITTER, 3, VCFHeaderLineType.Integer, READ_CONTEXT_JITTER_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_MAP_QUALITY, 2, VCFHeaderLineType.Integer, AVG_MAP_QUALITY_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_SUPPORT_DEPTH, 2, VCFHeaderLineType.Integer, "Raw allelic depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_SUPPORT_BASE_QUALITY, 2, VCFHeaderLineType.Integer, "Raw allelic base quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_DEPTH, 1, VCFHeaderLineType.Integer, "Raw read depth"));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_COUNT, VariantReadSupport.values().length, VCFHeaderLineType.Integer, READ_CONTEXT_COUNT_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_QUALITY, VariantReadSupport.values().length, VCFHeaderLineType.Integer, READ_CONTEXT_QUALITY_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_IMPROPER_PAIR, 1, VCFHeaderLineType.Integer, READ_CONTEXT_IMPROPER_PAIR_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(FRAG_STRAND_BIAS, 2, VCFHeaderLineType.Float, FRAG_STRAND_BIAS_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_STRAND_BIAS, 2, VCFHeaderLineType.Float, READ_STRAND_BIAS_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_BASE_QUAL, 1, VCFHeaderLineType.Integer, AVG_BASE_QUAL_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(TOTAL_RAW_BASE_QUAL, 1, VCFHeaderLineType.Integer, TOTAL_RAW_BASE_QUAL_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                UMI_TYPE_COUNTS, UMI_TYPE_COUNT, VCFHeaderLineType.Integer, UMI_TYPE_COUNTS_DESC));

        // filter options
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MNV_FILTER, "Filter duplicate MNV"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_SNV_MNV_FILTER, "Variant duplicate MNV vs SNV"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_INDEL_FILTER, "Variant duplicate SNV/MNV vs INDEL"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MIXED_GERMLINE_SOMATIC_FILTER, "Variant duplicate mixed somatic/germline"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MATCH, "Variant duplicate with different read contexts"));

        for(SoftFilter filter : SoftFilter.values())
        {
            header.addMetaDataLine(new VCFFilterHeaderLine(filter.filterName(), filter.vcfDescription()));
        }

        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        return header;
    }

    public static void appendHeader(final VCFHeader header)
    {
        if(!header.hasFormatLine(AVG_BASE_QUAL))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(AVG_BASE_QUAL, 1, VCFHeaderLineType.Integer, AVG_BASE_QUAL_DESC));
        }

        if(!header.hasFormatLine(AVG_MAP_QUALITY))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(AVG_MAP_QUALITY, 2, VCFHeaderLineType.Integer, AVG_MAP_QUALITY_DESC));
        }

        if(!header.hasFormatLine(UMI_TYPE_COUNTS))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(
                    UMI_TYPE_COUNTS, UMI_TYPE_COUNT, VCFHeaderLineType.Integer, UMI_TYPE_COUNTS_DESC));
        }

        if(!header.hasFormatLine(READ_STRAND_BIAS))
        {
            header.addMetaDataLine(new VCFFormatHeaderLine(READ_STRAND_BIAS, 1, VCFHeaderLineType.Float, READ_STRAND_BIAS_DESC));
        }

        if(!header.hasInfoLine(MAX_READ_EDGE_DISTANCE))
        {
            header.addMetaDataLine(new VCFInfoHeaderLine(MAX_READ_EDGE_DISTANCE, 1, VCFHeaderLineType.Integer, MAX_READ_EDGE_DISTANCE_DESC));
        }
    }

    @Override
    public void close()
    {
        mWriter.close();
    }
}
