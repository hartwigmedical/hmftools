package com.hartwig.hmftools.purple.somatic;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentPurple;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.RefGenomeData;
import com.hartwig.hmftools.purple.config.SomaticConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVCF {

    private final SomaticConfig somaticConfig;
    private final CommonConfig commonConfig;
    private final RefGenomeData refGenomeData;
    private final String inputVCF;
    private final String outputVCF;
    private boolean enabled;

    public SomaticVCF(final CommonConfig commonConfig, final SomaticConfig somaticConfig, final RefGenomeData refGenomeData) {
        this.commonConfig = commonConfig;
        this.outputVCF = commonConfig.outputDirectory() + File.separator + commonConfig.tumorSample() + ".purple.somatic.vcf.gz";
        this.enabled = somaticConfig.file().isPresent();
        this.inputVCF = enabled ? somaticConfig.file().get().toString() : "";
        this.refGenomeData = refGenomeData;
        this.somaticConfig = somaticConfig;
    }

    public void write(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions, @NotNull final List<PeakModel> somaticPeaks) throws IOException {

        if (enabled) {

            try (IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(refGenomeData.refGenome()));
                    VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false);
                    VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                            .setReferenceDictionary(vcfReader.getFileHeader().getSequenceDictionary())
                            .setIndexCreator(new TabixIndexCreator(vcfReader.getFileHeader().getSequenceDictionary(), new TabixFormat()))
                            .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                            .build()) {

                final VariantContextEnrichmentPurple enricher = new VariantContextEnrichmentPurple(somaticConfig.clonalityMaxPloidy(),
                        somaticConfig.clonalityBinWidth(),
                        commonConfig.version(),
                        commonConfig.tumorSample(),
                        indexedFastaSequenceFile,
                        purityAdjuster,
                        copyNumbers,
                        fittedRegions,
                        somaticPeaks,
                        writer::add);

                final VCFHeader header = enricher.enrichHeader(vcfReader.getFileHeader());
                writer.writeHeader(header);

                for (VariantContext context : vcfReader) {
                    enricher.accept(context);
                }

                enricher.flush();

            }
        }
    }
}
