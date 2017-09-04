package com.hartwig.hmftools.purple.baf;

import java.io.IOException;
import java.util.function.Supplier;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.baf.TumorBAF;
import com.hartwig.hmftools.common.baf.TumorBAFFile;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.purple.config.BAFConfig;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BAFSupplier implements Supplier<Multimap<String, TumorBAF>> {

    private static final Logger LOGGER = LogManager.getLogger(BAFSupplier.class);

    @NotNull
    private final Supplier<Multimap<String, TumorBAF>> supplier;

    public BAFSupplier(@NotNull final CommonConfig config, @NotNull final BAFConfig bafConfig)
            throws ParseException, IOException, HartwigException {

        if (bafConfig.bafFile().isPresent()) {
            final String bafFile = bafConfig.bafFile().get().toString();
            LOGGER.info("Reading BAFs from {}", bafFile);
            supplier = new FileSupplier(bafFile);
        } else {
            final String bafFile = TumorBAFFile.generatePurpleFilename(config.outputDirectory(), config.tumorSample());
            LOGGER.info("Generating BAFs from germline VCF");
            supplier = new VCFSupplier(bafConfig.bafVCFFile().get());

            LOGGER.info("Persisting BAFs to {}", bafFile);
            TumorBAFFile.write(bafFile, supplier.get());
        }
    }

    @Override
    @NotNull
    public Multimap<String, TumorBAF> get() {
        return supplier.get();
    }
}
