package com.hartwig.hmftools.sage.vcf;

import static com.hartwig.hmftools.common.variant.VariantHeader.PASS;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
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
    public static final String VERSION_META_DATA = "sageVersion";
    public static final String READ_CONTEXT = "RC";
    public static final String READ_CONTEXT_LEFT_FLANK = "RC_LF";
    public static final String READ_CONTEXT_RIGHT_FLANK = "RC_RF";
    public static final String READ_CONTEXT_INDEX = "RC_IDX";

    public static final String DEDUP_MNV_FILTER = "dedupMnv";
    public static final String DEDUP_MIXED_GERMLINE_SOMATIC_FILTER = "dedupMixedGermlineSomatic";
    public static final String DEDUP_SNV_MNV_FILTER = "dedupSnvMnv";
    public static final String DEDUP_INDEL_FILTER = "dedupIndel";
    public static final String DEDUP_MATCH = "dedupMatch";

    public static final String READ_CONTEXT_JITTER = "RC_JIT";
    private static final String READ_CONTEXT_JITTER_DESCRIPTION = "Read context jitter [Shortened, Lengthened, QualityPenalty]";

    public static final String READ_CONTEXT_EVENTS = "RC_NM";
    private static final String READ_CONTEXT_EVENTS_DESCRIPTION = "Minimum number of events in read";

    public static final String READ_CONTEXT_COUNT = "RC_CNT";
    private static final String READ_CONTEXT_COUNT_DESCRIPTION =
            "Read context counts [Full, Partial, Core, Realigned, Alt, Reference, Total]";
    public static final String READ_CONTEXT_REPEAT_COUNT = "RC_REPC";
    private static final String READ_CONTEXT_REPEAT_COUNT_DESCRIPTION = "Repeat count at read context";
    public static final String READ_CONTEXT_REPEAT_SEQUENCE = "RC_REPS";
    private static final String READ_CONTEXT_REPEAT_SEQUENCE_DESCRIPTION = "Repeat sequence at read context";
    public static final String READ_CONTEXT_MICRO_HOMOLOGY = "RC_MH";
    private static final String READ_CONTEXT_MICRO_HOMOLOGY_DESCRIPTION = "Micro-homology at read context";
    public static final String READ_CONTEXT_QUALITY = "RC_QUAL";
    private static final String READ_CONTEXT_QUALITY_DESCRIPTION =
            "Read context quality [Full, Partial, Core, Realigned, Alt, Reference, Total]";
    private static final String READ_CONTEXT_AF_DESCRIPTION =
            "Allelic frequency calculated from read context counts as (Full + Partial + Core + Realigned + Alt) / Coverage";

    public static final String READ_CONTEXT_IMPROPER_PAIR = "RC_IPC";
    private static final String READ_CONTEXT_IMPROPER_PAIR_DESCRIPTION = "Read context improper pair count";

    public static final String RAW_DEPTH = "RDP";
    public static final String RAW_ALLELIC_DEPTH = "RAD";
    public static final String RAW_ALLELIC_BASE_QUALITY = "RABQ";

    public static final String STRAND_BIAS = "SB";
    public static final String STRAND_BIAS_DESC = "Strand bias - percentage of first-in-pair reads";

    public static final String AVG_BASE_QUAL = "ABQ";
    public static final String AVG_BASE_QUAL_DESC = "Average calculated base quality";

    public static final String MIXED_SOMATIC_GERMLINE = "MSG";
    public static final String MIXED_SOMATIC_GERMLINE_DESCRIPTION = "Mixed Somatic and Germline variants";

    public static final String LOCAL_PHASE_SET_READ_COUNT = "LPS_RC";
    private static final String LPS_READ_COUNT_DESCRIPTION = "Local Phase Set Read Count";

    private final VariantContextWriter mWriter;
    private final Consumer<VariantContext> mConsumer;

    public VariantVCF(final IndexedFastaSequenceFile reference, final SageConfig config)
    {
        final SAMSequenceDictionary sequenceDictionary = reference.getSequenceDictionary();

        mWriter = new VariantContextWriterBuilder().setOutputFile(config.OutputFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(sequenceDictionary)
                .build();
        SomaticRefContextEnrichment enrichment = new SomaticRefContextEnrichment(reference, mWriter::add);
        mConsumer = enrichment;

        final VCFHeader header = enrichment.enrichHeader(header(config));

        final SAMSequenceDictionary condensedDictionary = new SAMSequenceDictionary();
        for(SAMSequenceRecord sequence : sequenceDictionary.getSequences())
        {
            if(HumanChromosome.contains(sequence.getContig()) || MitochondrialChromosome.contains(sequence.getContig()))
            {
                condensedDictionary.addSequence(sequence);
            }
        }

        header.setSequenceDictionary(condensedDictionary);
        mWriter.writeHeader(header);
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
        mConsumer = mWriter::add;
        mWriter.writeHeader(newHeader);
    }

    public void write(final VariantContext context)
    {
        mConsumer.accept(context);
    }

    private static VCFHeader header(final SageConfig config)
    {
        final List<String> samples = Lists.newArrayList();
        samples.addAll(config.ReferenceIds);
        samples.addAll(config.TumorIds);
        return header(config.Version, samples);
    }

    private static VCFHeader header(final String version, final List<String> allSamples)
    {
        VCFHeader header = SageMetaData.addSageMetaData(new VCFHeader(Collections.emptySet(), allSamples));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                LOCAL_PHASE_SET_READ_COUNT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, LPS_READ_COUNT_DESCRIPTION));

        header.addMetaDataLine(new VCFHeaderLine(VERSION_META_DATA, version));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        header.addMetaDataLine(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));
        header.addMetaDataLine(new VCFFormatHeaderLine(
                VCFConstants.ALLELE_FREQUENCY_KEY, 1, VCFHeaderLineType.Float, READ_CONTEXT_AF_DESCRIPTION));

        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_JITTER, 3, VCFHeaderLineType.Integer, READ_CONTEXT_JITTER_DESCRIPTION));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_ALLELIC_DEPTH, 2, VCFHeaderLineType.Integer, "Raw allelic depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_ALLELIC_BASE_QUALITY, 2, VCFHeaderLineType.Integer, "Raw allelic base quality"));
        header.addMetaDataLine(new VCFFormatHeaderLine(RAW_DEPTH, 1, VCFHeaderLineType.Integer, "Raw read depth"));
        header.addMetaDataLine(new VCFFormatHeaderLine(READ_CONTEXT_COUNT, 7, VCFHeaderLineType.Integer, READ_CONTEXT_COUNT_DESCRIPTION));
        // header.addMetaDataLine(new VCFFormatHeaderLine(SC_INSERT_SUPPORT, 1, VCFHeaderLineType.Integer, SC_INSERT_SUPPORT_DESCRIPTION));
        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_IMPROPER_PAIR, 1, VCFHeaderLineType.Integer, READ_CONTEXT_IMPROPER_PAIR_DESCRIPTION));
        header.addMetaDataLine(new VCFFormatHeaderLine(
                READ_CONTEXT_QUALITY, 7, VCFHeaderLineType.Integer, READ_CONTEXT_QUALITY_DESCRIPTION));

        header.addMetaDataLine(new VCFFormatHeaderLine(STRAND_BIAS, 1, VCFHeaderLineType.Float, STRAND_BIAS_DESC));
        header.addMetaDataLine(new VCFFormatHeaderLine(AVG_BASE_QUAL, 1, VCFHeaderLineType.Integer, AVG_BASE_QUAL_DESC));

        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_EVENTS, 1, VCFHeaderLineType.Integer, READ_CONTEXT_EVENTS_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT, 1, VCFHeaderLineType.String, "Read context core"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_INDEX, 1, VCFHeaderLineType.Integer, "Read context index"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_LEFT_FLANK, 1, VCFHeaderLineType.String, "Read context left flank"));
        header.addMetaDataLine(new VCFInfoHeaderLine(READ_CONTEXT_RIGHT_FLANK, 1, VCFHeaderLineType.String, "Read context right flank"));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_COUNT,1, VCFHeaderLineType.Integer, READ_CONTEXT_REPEAT_COUNT_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_REPEAT_SEQUENCE, 1, VCFHeaderLineType.String, READ_CONTEXT_REPEAT_SEQUENCE_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                READ_CONTEXT_MICRO_HOMOLOGY, 1, VCFHeaderLineType.String, READ_CONTEXT_MICRO_HOMOLOGY_DESCRIPTION));
        header.addMetaDataLine(new VCFInfoHeaderLine(
                MIXED_SOMATIC_GERMLINE, 1, VCFHeaderLineType.Integer, MIXED_SOMATIC_GERMLINE_DESCRIPTION));

        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MNV_FILTER, "Filter duplicate MNV"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_SNV_MNV_FILTER, "Variant duplicate MNV vs SNV"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_INDEL_FILTER, "Variant duplicate SNV/MNV vs INDEL"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MIXED_GERMLINE_SOMATIC_FILTER, "Variant duplicate mixed somatic/germline"));
        header.addMetaDataLine(new VCFFilterHeaderLine(DEDUP_MATCH, "Variant duplicate with different read contexts"));

        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MIN_AVG_BASE_QUALITY.filterName(), "Variant average base quality below limit"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.STRAND_BIAS.filterName(), "Variant exceeds strand bias limit"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MIN_TUMOR_QUAL.filterName(), "Insufficient tumor quality"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MIN_TUMOR_VAF.filterName(), "Insufficient tumor VAF"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MIN_GERMLINE_DEPTH.filterName(), "Insufficient germline depth"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MAX_GERMLINE_VAF.filterName(), "Excess germline VAF"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MAX_GERMLINE_REL_RAW_BASE_QUAL.filterName(),
                "Excess germline relative quality"));
        header.addMetaDataLine(new VCFFilterHeaderLine(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.filterName(), "Excess germline alt support"));
        header.addMetaDataLine(new VCFFilterHeaderLine(PASS, "All filters passed"));

        return header;
    }

    @Override
    public void close()
    {
        mWriter.close();
    }

}
