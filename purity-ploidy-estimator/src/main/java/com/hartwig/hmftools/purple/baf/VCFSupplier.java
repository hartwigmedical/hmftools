package com.hartwig.hmftools.purple.baf;

import java.io.IOException;
import java.util.List;
import java.util.function.Supplier;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.common.purple.baf.TumorBAFFactory;
import com.hartwig.hmftools.common.purple.baf.TumorBAFFile;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFGermlineFile;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

class VCFSupplier implements Supplier<Multimap<String, TumorBAF>> {

    private static final Logger LOGGER = LogManager.getLogger(VCFSupplier.class);

    private static final double MIN_REF_ALLELE_FREQUENCY = 0.4;
    private static final double MAX_REF_ALLELE_FREQUENCY = 0.65;
    private static final int MIN_COMBINED_DEPTH = 10;
    private static final int MAX_COMBINED_DEPTH = 100;

    private final Multimap<String, TumorBAF> bafs;

    VCFSupplier(final CommonConfig config, final String vcfExtension) throws IOException, HartwigException {

        final TumorBAFFactory factory =
                new TumorBAFFactory(MIN_REF_ALLELE_FREQUENCY, MAX_REF_ALLELE_FREQUENCY, MIN_COMBINED_DEPTH, MAX_COMBINED_DEPTH);

        LOGGER.info("Calculating BAFs from germline VCF");
        final VCFGermlineFile vcfFile = VCFFileLoader.loadGermlineVCF(config.runDirectory(), vcfExtension);
        final List<GermlineVariant> variants = VariantFilter.passOnly(vcfFile.variants());
        bafs = factory.createBAF(variants);

        final String filename = TumorBAFFile.generateFilename(config.outputDirectory(), config.tumorSample());
        LOGGER.info("Persisting BAFs to {}", filename);
        TumorBAFFile.write(filename, bafs);
    }

    @Override
    public Multimap<String, TumorBAF> get() {
        return bafs;
    }
}
