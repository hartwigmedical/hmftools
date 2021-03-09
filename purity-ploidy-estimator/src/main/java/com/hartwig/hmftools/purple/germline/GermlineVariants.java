package com.hartwig.hmftools.purple.germline;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.VariantHeader;
import com.hartwig.hmftools.common.variant.enrich.GermlineVariantEnrichment;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.ConfigSupplier;
import com.hartwig.hmftools.purple.config.RefGenomeData;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVariants {
    private static final Logger LOGGER = LogManager.getLogger(GermlineVariants.class);

    private final ConfigSupplier configSupplier;
    private final CommonConfig commonConfig;
    private final RefGenomeData refGenomeData;
    private final String outputVCF;
    private final List<CanonicalTranscript> transcripts;
    private final DriverGenePanel genePanel;
    private final List<VariantContext> reportableVariants;

    public GermlineVariants(final ConfigSupplier configSupplier) {
        this.genePanel = configSupplier.driverCatalogConfig().genePanel();
        this.commonConfig = configSupplier.commonConfig();
        this.refGenomeData = configSupplier.refGenomeConfig();
        this.configSupplier = configSupplier;
        this.outputVCF = commonConfig.outputDirectory() + File.separator + commonConfig.tumorSample() + ".purple.germline.vcf.gz";
        this.transcripts =
                configSupplier.refGenomeConfig().isHg38() ? CanonicalTranscriptFactory.create38() : CanonicalTranscriptFactory.create37();
        this.reportableVariants = Lists.newArrayList();
    }

    @NotNull
    public List<VariantContext> reportableVariants() {
        return reportableVariants;
    }

    public void processAndWrite(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<VariantContext> somaticVariants) throws IOException {

        final Optional<File> optionalInputVCF = configSupplier.germlineConfig().file();
        if (optionalInputVCF.isPresent()) {

            LOGGER.info("Loading germline variants from {}", optionalInputVCF.get());
            LOGGER.info("Enriching germline variants");

            try (IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomeData.refGenome()));
                    VCFFileReader vcfReader = new VCFFileReader(optionalInputVCF.get(), false);
                    VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                            .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                            .build()) {

                final Consumer<VariantContext> consumer = context -> {
                    if (context.getAttributeAsBoolean(VariantHeader.REPORTED_FLAG, false)) {
                        reportableVariants.add(context);
                    }
                    writer.add(context);
                };

                final GermlineVariantEnrichment enrichment = new GermlineVariantEnrichment(commonConfig.version(),
                        commonConfig.refSample(),
                        commonConfig.tumorSample(),
                        indexedFastaSequenceFile,
                        purityAdjuster,
                        copyNumbers,
                        genePanel,
                        transcripts,
                        configSupplier.driverCatalogConfig().germlineHotspots(),
                        somaticVariants,
                        consumer);

                writer.writeHeader(enrichment.enrichHeader(vcfReader.getFileHeader()));
                for (VariantContext context : vcfReader) {
                    enrichment.accept(context);
                }
                enrichment.flush();
            }
        }
    }
}
