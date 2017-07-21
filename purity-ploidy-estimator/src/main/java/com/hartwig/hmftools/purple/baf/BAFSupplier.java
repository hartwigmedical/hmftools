package com.hartwig.hmftools.purple.baf;

import static com.hartwig.hmftools.purple.CommandLineUtil.defaultValue;

import java.io.IOException;
import java.util.function.Supplier;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.baf.TumorBAF;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class BAFSupplier implements Supplier<Multimap<String, TumorBAF>> {

    private static final String BAF_FILE = "baf_file"; //TODO: Add support
    private static final String VCF_EXTENSION = "vcf_extension";
    private static final String VCF_EXTENSION_DEFAULT = ".annotated_sliced.vcf";

    public static void addOptions(Options options) {
        options.addOption(VCF_EXTENSION, true, "VCF file extension. Defaults to " + VCF_EXTENSION_DEFAULT);
    }

    private final Supplier<Multimap<String, TumorBAF>> supplier;

    public BAFSupplier(final CommonConfig config, final CommandLine cmd) throws ParseException, IOException, HartwigException {
        final String vcfExtension = defaultValue(cmd, VCF_EXTENSION, VCF_EXTENSION_DEFAULT);
        supplier = new VCFBAFSupplier(config, vcfExtension);
    }

    @Override
    public Multimap<String, TumorBAF> get() {
        return supplier.get();
    }
}
