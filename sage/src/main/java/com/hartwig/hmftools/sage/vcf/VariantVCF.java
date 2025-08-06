package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_BASE_QUAL_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_EDGE_DISTANCE_PERC_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MICROHOMOLOGY_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MIN_COORDS_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RAW_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_RAW_BASE_QUAL_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_MICROHOMOLOGY_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.REPEAT_SEQUENCE;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TINC_RECOVERED_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TINC_RECOVERED_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TRINUCLEOTIDE_CONTEXT_DESC;
import static com.hartwig.hmftools.common.variant.SageVcfTags.CONSENSUS_TAG_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_BASE_QUAL_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_ALT_MAP_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MODIFIED_ALT_MAP_QUAL_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LOCAL_PHASE_SET_READ_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.LPS_READ_COUNT_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MIXED_SOMATIC_GERMLINE_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_AF_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_EVENTS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_STRAND_BIAS_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.SIMPLE_ALT_COUNT;
import static com.hartwig.hmftools.sage.vcf.VcfTags.SIMPLE_ALT_COUNT_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TUMOR_QUALITY_PROB;
import static com.hartwig.hmftools.sage.vcf.VcfTags.TUMOR_QUALITY_PROB_DESC;
import static com.hartwig.hmftools.sage.vcf.VcfTags.VERSION_META_DATA;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.common.variant.pon.PonCache;
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

    public VariantVCF(
            final IndexedFastaSequenceFile reference, final String sageVersion, final List<String> sampleIds, final String outputVcf,
            boolean runTinc)
    {
        final SAMSequenceDictionary sequenceDictionary = reference.getSequenceDictionary();

        mWriter = new VariantContextWriterBuilder().setOutputFile(outputVcf)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        VCFHeader vcfHeader = createHeader(sageVersion, sampleIds, runTinc);

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

    public VariantVCF(
            final IndexedFastaSequenceFile reference, final List<String> additionalSampleIds, final VCFHeader existingHeader,
            final String outputVcf)
    {
        Set<VCFHeaderLine> headerLines = existingHeader.getMetaDataInInputOrder();
        List<String> samples = Lists.newArrayList(existingHeader.getGenotypeSamples());
        samples.addAll(additionalSampleIds);

        VCFHeader newHeader = new VCFHeader(headerLines, samples);

        mWriter = new VariantContextWriterBuilder().setOutputFile(outputVcf)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(reference.getSequenceDictionary())
                .build();
        mWriter.writeHeader(newHeader);
    }

    public void write(final VariantContext context)
    {
        mWriter.add(context);
    }

    public static VCFHeader createHeader(final String version, final List<String> allSamples, final boolean runTinc)
    {
        VCFHeader header = new VCFHeader(Collections.emptySet(), allSamples);

        header.addMetaDataLine(new VCFInfoHeaderLine(TIER, 1, VCFHeaderLineType.String, TIER_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, LOCAL_PHASE_SET_DESC));

        // standard fields
        header.addMetaDataLine(new VCFHeaderLine(VERSION_META_DATA, version));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));

        // info fields
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_EVENTS, 1, VCFHeaderLineType.Integer, READ_CONTEXT_EVENTS_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_INFO, 1, VCFHeaderLineType.String, READ_CONTEXT_INFO_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(TRINUCLEOTIDE_CONTEXT, 1, VCFHeaderLineType.String, TRINUCLEOTIDE_CONTEXT_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_SEQUENCE, 1, VCFHeaderLineType.String, READ_CONTEXT_REPEAT_SEQUENCE_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_COUNT, 1, VCFHeaderLineType.Integer, READ_CONTEXT_REPEAT_COUNT_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_SEQUENCE, 1, VCFHeaderLineType.String, REPEAT_SEQUENCE_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(REPEAT_COUNT, 1, VCFHeaderLineType.Integer, REPEAT_COUNT_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(MICROHOMOLOGY, 1, VCFHeaderLineType.String, MICROHOMOLOGY_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_MICROHOMOLOGY, 1, VCFHeaderLineType.String, READ_CONTEXT_MICROHOMOLOGY_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                MIXED_SOMATIC_GERMLINE, 1, VCFHeaderLineType.Integer, MIXED_SOMATIC_GERMLINE_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET_READ_COUNT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, LPS_READ_COUNT_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(MAX_READ_EDGE_DISTANCE, 1, VCFHeaderLineType.Float, MAX_READ_EDGE_DISTANCE_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(NEARBY_INDEL_FLAG, 0, VCFHeaderLineType.Flag, NEARBY_INDEL_FLAG_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(TUMOR_QUALITY_PROB, 1, VCFHeaderLineType.Float, TUMOR_QUALITY_PROB_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(MAP_QUAL_FACTOR, 1, VCFHeaderLineType.Float, MAP_QUAL_FACTOR_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(TINC_RECOVERED_FLAG, 0, VCFHeaderLineType.Flag, TINC_RECOVERED_DESC));

        // genotype fields
        addGenotypeHeader(header);

        for(SoftFilter filter : SoftFilter.values())
        {
            header.addMetaDataLine(new VCFFilterHeaderLine(filter.filterName(), filter.vcfDescription()));
        }

        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        if(runTinc)
        {
            GnomadCache.addAnnotationHeader(header);
            PonCache.addAnnotationHeader(header);
        }

        return header;
    }

    public static void addGenotypeHeader(final VCFHeader header)
    {
        // call from Sage append as well for new samples, and this handles the additional in later versions of new genotype fields

        // add in alphabetical order
        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_BASE_QUAL, 2, VCFHeaderLineType.Integer, AVG_BASE_QUAL_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_EDGE_DISTANCE_PERC, 2, VCFHeaderLineType.Float, AVG_EDGE_DISTANCE_PERC_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                VCFConstants.ALLELE_FREQUENCY_KEY, 1, VCFHeaderLineType.Float, READ_CONTEXT_AF_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                AVG_MODIFIED_BASE_QUAL, 1, VCFHeaderLineType.Integer, AVG_MODIFIED_BASE_QUAL_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                AVG_MODIFIED_ALT_MAP_QUAL, 1, VCFHeaderLineType.Integer, AVG_MODIFIED_ALT_MAP_QUAL_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_MAP_QUALITY, 2, VCFHeaderLineType.Integer, AVG_MAP_QUALITY_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(MIN_COORDS_COUNT, 1, VCFHeaderLineType.Integer, MIN_COORDS_COUNT_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_RAW_BASE_QUAL, 1, VCFHeaderLineType.Integer, AVG_RAW_BASE_QUAL_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_COUNT, VariantReadSupport.values().length, VCFHeaderLineType.Integer, READ_CONTEXT_COUNT_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_IMPROPER_PAIR, 1, VCFHeaderLineType.Integer, READ_CONTEXT_IMPROPER_PAIR_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_JITTER, 3, VCFHeaderLineType.Integer, READ_CONTEXT_JITTER_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_QUALITY, VariantReadSupport.values().length, VCFHeaderLineType.Integer, READ_CONTEXT_QUALITY_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(READ_STRAND_BIAS, 2, VCFHeaderLineType.Float, READ_STRAND_BIAS_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(SIMPLE_ALT_COUNT, 1, VCFHeaderLineType.Integer, SIMPLE_ALT_COUNT_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(FRAG_STRAND_BIAS, 2, VCFHeaderLineType.Float, FRAG_STRAND_BIAS_DESC));

        header.addMetaDataLine(new VCFFormatHeaderLine(UMI_TYPE_COUNTS, CONSENSUS_TAG_TYPE_COUNT, VCFHeaderLineType.Integer, UMI_TYPE_COUNTS_DESC));
    }

    @Override
    public void close()
    {
        mWriter.close();
    }
}
