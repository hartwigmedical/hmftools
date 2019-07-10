package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_GERMLINE_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_MINOR_ALLELE_PLOIDY_INFO;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PURPLE_PLOIDY_INFO;

import java.io.File;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
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
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SomaticVCF {

    private static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    private static final String PURPLE_AF_DESC = "Purity adjusted allelic frequency of variant";
    private static final String PURPLE_PLOIDY_DESC = "Purity adjusted ploidy of variant";
    private static final String PURPLE_BIALLELIC_DESC = "Variant is biallelic";

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

            final PurityAdjustedSomaticVariantFactory enricher =
                    new PurityAdjustedSomaticVariantFactory(commonConfig.tumorSample(), purityAdjuster, copyNumbers, fittedRegions);

            final VCFFileReader vcfReader = new VCFFileReader(new File(inputVCF), false);
            final VCFHeader header = generateOutputHeader(commonConfig.version(), vcfReader.getFileHeader());

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

    @NotNull
    @VisibleForTesting
    private static VCFHeader generateOutputHeader(@NotNull final String purpleVersion, @NotNull final VCFHeader template) {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_PLOIDY_INFO, 1, VCFHeaderLineType.Float, PURPLE_PLOIDY_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_PLOIDY_INFO,
                1,
                VCFHeaderLineType.Float,
                PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(SomaticVariantFactory.PURPLE_BIALLELIC_FLAG,
                0,
                VCFHeaderLineType.Flag,
                PURPLE_BIALLELIC_DESC));


        return outputVCFHeader;
    }

}
