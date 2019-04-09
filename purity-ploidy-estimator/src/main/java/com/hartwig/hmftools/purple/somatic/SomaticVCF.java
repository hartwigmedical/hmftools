package com.hartwig.hmftools.purple.somatic;

import java.io.File;
import java.util.List;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.purple.config.CommonConfig;
import com.hartwig.hmftools.purple.config.SomaticConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVCF {

    private final CommonConfig commonConfig;
    private final String inputVCF;
    private final String outputVCF;
    private boolean enabled;

    public SomaticVCF(final CommonConfig commonConfig, final SomaticConfig somaticConfig) {
        this.commonConfig = commonConfig;
        this.outputVCF = commonConfig.outputDirectory() + File.separator + commonConfig.tumorSample() + ".purple.somatic.vcf.gz";
        this.enabled = somaticConfig.file().isPresent();
        this.inputVCF = enabled ? somaticConfig.file().get().toString() : "";

    }

    public void write(@NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions) {

        if (enabled) {

            final SomaticEnrichment enricher =
                    new SomaticEnrichment(purityAdjuster, copyNumbers, fittedRegions, commonConfig.tumorSample());
            final VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false);
            final VCFHeader header = SomaticEnrichment.generateOutputHeader(commonConfig.version(), vcfReader.getFileHeader());

            final VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(outputVCF)
                    .setReferenceDictionary(header.getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(header.getSequenceDictionary(), new TabixFormat()))
                    .setOption(htsjdk.variant.variantcontext.writer.Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();
            writer.writeHeader(header);

            for (VariantContext context : vcfReader) {
                writer.add(enricher.enrich(context));
            }

            vcfReader.close();
            writer.close();
        }
    }

}
