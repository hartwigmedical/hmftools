package com.hartwig.hmftools.serve.extraction.util;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public final class VCFWriterFactory {

    public static final String INPUT_FIELD = "input";
    public static final String SOURCES_FIELD = "sources";

    private VCFWriterFactory() {
    }

    @NotNull
    public static VariantContextWriter openIndexedVCFWriter(@NotNull String outputVcf, @NotNull IndexedFastaSequenceFile refSequence,
            @NotNull String sources) {
        SAMSequenceDictionary sequenceDictionary = refSequence.getSequenceDictionary();
        VariantContextWriter writer = createBaseWriterBuilder(outputVcf).modifyOption(Options.INDEX_ON_THE_FLY, true)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        VCFHeader header = createBaseHeader(sources);

        SAMSequenceDictionary condensedDictionary = new SAMSequenceDictionary();
        for (SAMSequenceRecord sequence : sequenceDictionary.getSequences()) {
            if (HumanChromosome.contains(sequence.getContig()) || MitochondrialChromosome.contains(sequence.getContig())) {
                condensedDictionary.addSequence(sequence);
            }
        }

        header.setSequenceDictionary(condensedDictionary);
        writer.writeHeader(header);
        return writer;
    }

    @NotNull
    public static VariantContextWriter openVCFWriter(@NotNull String outputVcf, @NotNull String sources) {
        VariantContextWriter writer = createBaseWriterBuilder(outputVcf).modifyOption(Options.INDEX_ON_THE_FLY, false).build();
        writer.writeHeader(createBaseHeader(sources));
        return writer;
    }

    @NotNull
    private static VariantContextWriterBuilder createBaseWriterBuilder(@NotNull String outputVcf) {
        return new VariantContextWriterBuilder().setOutputFile(outputVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .modifyOption(Options.USE_ASYNC_IO, false);
    }

    @NotNull
    private static VCFHeader createBaseHeader(@NotNull String sources) {
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList());
        header.addMetaDataLine(new VCFInfoHeaderLine(INPUT_FIELD,
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "Input based on which SERVE generated this record"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SOURCES_FIELD,
                VCFHeaderLineCount.UNBOUNDED,
                VCFHeaderLineType.String,
                "SERVE sources that contained this record [" + sources + "]"));
        return header;
    }
}
