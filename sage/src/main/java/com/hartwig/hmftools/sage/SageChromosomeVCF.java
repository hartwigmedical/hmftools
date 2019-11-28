package com.hartwig.hmftools.sage;

import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.sage.config.SageConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SageChromosomeVCF implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SageChromosomeVCF.class);

    private final String filename;
    private final VariantContextWriter writer;

    public SageChromosomeVCF(@NotNull final String chromosome, @NotNull final SageConfig config) throws IOException {
        filename = File.createTempFile("sage." + chromosome + ".", ".vcf.gz").toString();
        writer = new VariantContextWriterBuilder().setOutputFile(filename)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .modifyOption(Options.USE_ASYNC_IO, true)
                .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, true)
                .build();
        final VCFHeader header = SageVCF.header(config.reference(), config.tumor());

        LOGGER.info("Creating temporary file: {}", filename);
        writer.writeHeader(header);
    }

    @NotNull
    public String filename() {
        return filename;
    }

    public void write(@NotNull final VariantContext context) {
        writer.add(context);
    }

    @Override
    public void close() {
        writer.close();
    }

}
