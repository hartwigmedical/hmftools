package com.hartwig.hmftools.purple.somatic;

import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.purple.PurityPloidyEstimateApplication;
import com.hartwig.hmftools.purple.config.SomaticConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SomaticSupplier implements Supplier<List<SomaticVariant>> {

    private static final Logger LOGGER = LogManager.getLogger(PurityPloidyEstimateApplication.class);

    private final SomaticConfig config;

    public SomaticSupplier(final SomaticConfig config) {
        this.config = config;
    }

    @Override
    public List<SomaticVariant> get() {
        if (config.file().isPresent()) {
            String filename = config.file().toString();
            try {
                List<SomaticVariant> result = VCFFileLoader.loadSomaticVCF(filename)
                        .variants()
                        .stream()
                        .filter(x -> x.type() == VariantType.SNP)
                        .filter(VariantFilter::isPass)
                        .collect(Collectors.toList());
                LOGGER.info("Reading somatic variants from {}", filename);
                return result;
            } catch (Exception e) {
                return Collections.emptyList();
            }
        }

        return Collections.emptyList();
    }

}

