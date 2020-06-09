package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PurityEnrichment implements VariantContextEnrichment {

    public static final String PURPLE_CN_INFO = "PURPLE_CN";

    @Deprecated
    public static final String PURPLE_MINOR_ALLELE_PLOIDY_INFO = "PURPLE_MAP";
    public static final String PURPLE_MINOR_ALLELE_CN_INFO = "PURPLE_MACN";

    @Deprecated
    public static final String PURPLE_VARIANT_PLOIDY_INFO = "PURPLE_PLOIDY";
    public static final String PURPLE_VARIANT_CN_INFO = "PURPLE_VCN";

    public static final String PURPLE_AF_INFO = "PURPLE_AF";
    public static final String PURPLE_GERMLINE_INFO = "PURPLE_GERMLINE";
    public static final String PURPLE_BIALLELIC_FLAG = "BIALLELIC";

    private static final String PURPLE_CN_DESC = "Purity adjusted copy number surrounding variant location";
    private static final String PURPLE_MINOR_ALLELE_PLOIDY_DESC = "Purity adjusted minor allele ploidy surrounding variant location";
    private static final String PURPLE_GERMLINE_DESC = "Germline classification surrounding variant location";

    private static final String PURPLE_AF_DESC = "Purity adjusted variant allelic frequency";
    private static final String PURPLE_PLOIDY_DESC = "Purity adjusted variant copy number";
    private static final String PURPLE_BIALLELIC_DESC = "Variant is biallelic";

    private final String purpleVersion;
    private final PurityAdjustedSomaticVariantFactory factory;
    private final Consumer<VariantContext> consumer;

    PurityEnrichment(@NotNull final String purpleVersion, @NotNull final String tumorSample,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final List<FittedRegion> fittedRegions, @NotNull final Consumer<VariantContext> consumer) {
        this.purpleVersion = purpleVersion;
        this.consumer = consumer;
        this.factory = new PurityAdjustedSomaticVariantFactory(tumorSample, purityAdjuster, copyNumbers, fittedRegions);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        consumer.accept(factory.enrich(context));
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF_INFO, 1, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_VARIANT_CN_INFO, 1, VCFHeaderLineType.Float, PURPLE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_MINOR_ALLELE_CN_INFO, 1,  VCFHeaderLineType.Float, PURPLE_MINOR_ALLELE_PLOIDY_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_GERMLINE_INFO, 1, VCFHeaderLineType.String, PURPLE_GERMLINE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_BIALLELIC_FLAG, 0, VCFHeaderLineType.Flag, PURPLE_BIALLELIC_DESC));

        return template;
    }

    @Override
    public void flush() {
        // Empty
    }
}
